
import os, os.path, sys
import glob
import hgsc_vcf
import logging
from collections import *
from collections import defaultdict, OrderedDict

logger = logging.getLogger()
logger.setLevel(logging.DEBUG)
ch = logging.StreamHandler()
ch.setLevel(logging.ERROR)
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
ch.setFormatter(formatter)
logger.addHandler(ch)

fh = logging.FileHandler("merge_debug.log")
fh.setLevel(logging.DEBUG)
fh.setFormatter(formatter)
logger.addHandler(fh)

##
# Goal:
#   Merge multiple VCF files from independent callers into a single unified file.
#
# Components required for the merging process:
#   1. Reader classes for chunked iteration over multiple VCFs
#   2. Record-merging logic to resolve overlapping variants
#   3. Writer (hgsc_vcf.Writer) to generate the combined VCF output
#   4. Header harmonization mechanism to unify metadata across callers

class MetaRecord(object):
    def __init__(self, caller, record):
        self.caller = caller
        self.record = record
    
    def _seq_pos(self, record):
        # Remove the prefix
        chrom = record['CHROM']
        if chrom.startswith("chr"):
            chrom = chrom[3:]

        try:

            if chrom == 'X':
                return 30 * 1000000000 + record['POS']
            elif chrom == 'Y':
                return 40 * 1000000000 + record['POS']
            elif chrom in ('M', 'MT'):
                return 50 * 1000000000 + record['POS']
            else:
                return int(chrom) * 1000000000 + record['POS']

        except ValueError:
            # For non-standard contigs
            return hash(record['CHROM']) % (10**12) + record['POS']

    def __cmp__(self, other):
        si = self._seq_pos(self.record)
        oi = self._seq_pos(other.record)
        return si - oi
        '''
        block commented out because we don't actually care that the alts aren't the same???
        if si != oi:
            return si - oi
        else:
            return sum([hash(a) for a in self.record['ALT']]) - sum([hash(a) for a in other.record['ALT']])
        '''
    def __repr__(self):
        return "<record chr=%s, pos=%s, alt=%s>" % (self.record['CHROM'], self.record['POS'], self.record['ALT'])

class MetaReader(object):
    def __init__(self, fobj):
        self.reader = hgsc_vcf.Reader(fobj)
        self.caller = fobj.name
        # get the normal and primary sample ids
        sampleMapping = {l.fields.get('ID'):l.fields.get('SampleTCGABarcode') for l in self.reader.header.get_headers('SAMPLE')}
        if 'PRIMARY' not in sampleMapping and 'METASTATIC' in sampleMapping:
            sampleMapping['PRIMARY'] = sampleMapping['METASTATIC']
        elif 'PRIMARY' not in sampleMapping and 'RECURRANCE' in sampleMapping:
            sampleMapping['PRIMARY'] = sampleMapping['RECURRANCE']
        logger.info("Sample mapping for %s: %s", fobj.name, sampleMapping)
        self.normal = sampleMapping.get('NORMAL', "Unknown normal sample ID")
        self.primary = sampleMapping.get('PRIMARY', "Unknown primary sample ID")
        self._next = None
        self.take() # call to take this time will return None but will also fast forward the reader to the next position
    def __cmp__(self, other):
        return self._next.__cmp__(other._next)

    def peek(self):
        return self._next
    
    def take(self):
        old = self._next
        # fast forward to the next
        new = None
        while True:
            try:
                n = self.reader.next()
#                if 'GL' in n['CHROM']:
#                    logger.info("GL in chrom %s", n['CHROM'])
#                    logger.info("Closing %s", self.caller)
                    #self.reader.fobj.close()
#                    new = None
#                else:
                 #   if 'NORMAL' not in n['SAMPLES']:
                 #       n['SAMPLES']['NORMAL'] = n['SAMPLES'][self.normal]

                if 'PRIMARY' not in n['SAMPLES']:
                    if self.primary in n['SAMPLES']:
                        n['SAMPLES']['PRIMARY'] = n['SAMPLES'][self.primary]
                    elif 'METASTATIC' in n['SAMPLES']:
                        n['SAMPLES']['PRIMARY'] = n['SAMPLES']['METASTATIC']
                    elif 'RECURRANCE' in n['SAMPLES']:
                        n['SAMPLES']['PRIMARY'] = n['SAMPLES']['RECURRANCE']
                    else:
                        raise ValueError("Can't find the PRIMARY sample")

                #   if n['SAMPLES']['NORMAL']['GT'][0] not in  ('0/0', '0', '.', './.'):
                #        continue

                if 'PASS' not in n['FILTER']:
                    continue

                new = n
                break

            except StopIteration: # swallow the error and just set to None
                logger.info("Stopped iteration")
                logger.info("Closing %s", self.caller)
                self.reader.fobj.close()
                new = None
                break
        if new is None:
            self._next = None
        else:
            self._next = MetaRecord(self.caller, new)
        return old

    def __repr__(self):
        if self._next is None:
            return "<%s, closed>" % self.caller
        else:
            return "<%s, chr=%s, pos=%s>" % (self.caller, self._next.record['CHROM'], self._next.record['POS'])

    

class MultiVCFReader(object):
    def __init__(self, infiles, outfile, keys):
        self.buffer = 10
        self.infiles = {f:MetaReader(open(f, 'r')) for f in infiles}
        self.outfile = outfile
        self.outwriter = hgsc_vcf.Writer(open(self.outfile, 'w'), self.generate_header())

        # get the normal and primary sample ids
        sampleMapping = {l.fields.get('ID'):l.fields.get('SampleTCGABarcode') for l in self.infiles.values()[0].reader.header.get_headers('SAMPLE')}

        if 'PRIMARY' not in sampleMapping and 'METASTATIC' in sampleMapping:
            sampleMapping['PRIMARY'] = sampleMapping['METASTATIC']
        elif 'PRIMARY' not in sampleMapping and 'RECURRANCE' in sampleMapping:
            sampleMapping['PRIMARY'] = sampleMapping['RECURRANCE']

        self.normal = sampleMapping.get('NORMAL', "Unknown normal sample ID")
        self.primary = sampleMapping.get('PRIMARY', "Unknown primary sample ID")
        self.keymap = dict(zip(infiles, keys))

    # lets make this a generator so that we can keep up with the sorting
    # very likely that once something is in sorted order we can keep it that way
    def get_next(self):
        while True:
            nextsort = sorted([mr for mr in self.infiles.values() if mr.peek() is not None])
            if len(nextsort) < 1:
                raise StopIteration() # we can break here
            _raises = True
            for r in self._get_sorted_next(nextsort):
                _raises = False
                yield r
            if _raises:
                logger.info("Sorted %s", nextsort)
                if len(nextsort) > 1:
                    logger.info("n0 <= n1: %s", nextsort[0] <= nextsort[1])
                    logger.info("nextsort[0]._next: %s", nextsort[0]._next)
                raise ValueError("Entered sorted but did not yield")
    # helper generator 
    def _get_sorted_next(self, nextsort):
        if len(nextsort) < 2:
            n0 = nextsort[0]
            while n0.peek() is not None:
                yield n0.take()
        else:
            n0 = nextsort[0]
            n1 = nextsort[1]
            while n0.peek() is not None and n0 <= n1:
                yield n0.take()
            if n0.peek() is None:
                logger.info("Reached the end of %s", n0.caller)


    def generate_header(self):
        newHeader = hgsc_vcf.VCFHeader()
        newHeader.samples = ['PRIMARY'] # deterministic sample names now
        for infile in self.infiles.values():
            newHeader.headers += infile.reader.header.headers # append all of the headers together, who cares, we can sort out later
        newHeader.add_header('##COMMAND=<ID=vcf-merge>')
        return newHeader

    ##
    # yields batches of MetaRecord's
    def chunk(self):
        batch = []
        cpos = 0
        for r in self.get_next(): # run through the generator
            if len(batch) < 1:
                batch.append(r)
                cpos = r.record['POS'] + max(len(r.record['REF']), max([len(a) for a in r.record['ALT']])) - 1 + self.buffer
            elif r.record['CHROM'] != batch[0].record['CHROM']:
                yield batch
                batch = [r]
            elif r.record['POS'] > (batch[-1].record['POS'] + self.buffer):
                yield batch
                batch = [r]
                cpos = r.record['POS'] + max(len(r.record['REF']), max([len(a) for a in r.record['ALT']])) - 1 + self.buffer
            else:
                batch.append(r)
                cpos = max(cpos, r.record['POS'] + max(len(r.record['REF']), max([len(a) for a in r.record['ALT']])) - 1 + self.buffer)
        yield batch # yield the last batch

    # make this an iterable
    def __iter__(self):
        return self.chunk()


def parseInfo(merge, sample_type):

    # Merges FORMAT fields (GT, DP, AD, AF) for overlapping variants

    total_ref = total_alt = 0.0
    valid_count = 0
    gt_set = set()

    for m in merge:
        s = m.record['SAMPLES'][sample_type]
        raw_gt = s.get('GT', ['.'])[0]
        gt_set.add('0/1' if m.caller == 'Mutect2' else raw_gt)

        try:
            ad = s.get('AD', [])
            if len(ad) >= 2:
                ref = int(float(ad[0]))
                alt = int(float(ad[1]))
            else:
                continue

            total_ref += ref
            total_alt += alt
            valid_count += 1

        except Exception as e:
            raise RuntimeError("Error parsing sample info: {}".format(e))

    if valid_count == 0:
        return {'GT': ['.'], 'DP': ['0'], 'AD': ['0', '0'], 'AF': ['0.0000']}

    mean_ref = int(round(total_ref / valid_count))
    mean_alt = int(round(total_alt / valid_count))
    dp = mean_ref + mean_alt
    af = "%.4f" % (float(mean_alt) / dp) if dp > 0 else '0.0000'

    if len(gt_set) == 1:
        final_gt = list(gt_set)[0]
    elif '0/1' in gt_set:
        final_gt = '0/1'
    else:
        final_gt = list(gt_set)[0]

    return {
        'GT': [final_gt],
        'DP': [str(dp)],
        'AD': [str(mean_ref), str(mean_alt)],
        'AF': [af]
    }

def resolve_merge(merge, callermap):
    results = []
    grouped = defaultdict(list)

    for m in merge:
        alt = str(m.record['ALT'][0])
        key = (m.record['CHROM'], m.record['POS'], m.record['REF'], alt)
        grouped[key].append(m)

    for (chrom, pos, ref, alt), group in grouped.items():
        ref = str(group[0].record['REF'])
        unique_callers = set(m.caller for m in group)
        # Non-overlapping variant: only one caller reports this site
        if len(unique_callers) == 1:
            orig = group[0].record.copy()
            s = orig['SAMPLES']['PRIMARY']
            ad = s.get('AD', ['0', '0'])
            if len(ad) >= 2:
                ref_count = int(float(ad[0]))
                alt_count = int(float(ad[1]))
            else:
                ref_count = alt_count = 0

            dp_int = ref_count + alt_count
            af = "%.4f" % (float(alt_count) / dp_int) if dp_int > 0 else '0.0000'
            ad_out = [str(ref_count), str(alt_count)]
            dp = str(dp_int)

            caller = group[0].caller
            raw_gt = s.get('GT', ['.'])[0]
            gt = '0/1' if caller == 'Mutect2' else raw_gt

            new_record = OrderedDict()
            new_record['CHROM'] = chrom
            new_record['POS'] = pos
            new_record['ID'] = '.'
            new_record['REF'] = ref
            new_record['ALT'] = [alt]
            new_record['QUAL'] = str(orig.get('QUAL', '.'))
            new_record['FILTER'] = [str(f) for f in orig.get('FILTER', ['PASS'])]
            new_record['INFO'] = {'CALLERS': callermap[caller]}
            new_record['FORMAT'] = ['GT', 'DP', 'AD', 'AF']
            new_record['SAMPLES'] = OrderedDict()
            new_record['SAMPLES']['PRIMARY'] = {
                'GT': [gt],
                'DP': [dp],
                'AD': ad_out,
                'AF': [af]
            }

            results.append((new_record, callermap[caller]))

        # Overlapping variant: multiple callers report this site, merge sample fields
        else:
            merged_info = parseInfo(group, 'PRIMARY')

            new_record = OrderedDict()
            new_record['CHROM'] = chrom
            new_record['POS'] = pos
            new_record['ID'] = '.'
            new_record['REF'] = ref
            new_record['ALT'] = [alt]
            vardict_records = [m for m in group if callermap[m.caller] == "VarDict"]
            new_record['QUAL'] = vardict_records[0].record.get('QUAL', '.') if vardict_records else group[0].record.get('QUAL', '.')
            new_record['FILTER'] = ['PASS']
            new_record['INFO'] = {'CALLERS': '|'.join(sorted(callermap[m.caller] for m in group))}
            new_record['FORMAT'] = ['GT', 'DP', 'AD', 'AF']
            new_record['SAMPLES'] = OrderedDict()
            new_record['SAMPLES']['PRIMARY'] = merged_info

            results.append((new_record, new_record['INFO']['CALLERS']))

    return results



def contains_pindel(batch):
    for r in batch:
        if 'pindel' in r.caller:
            return True
    return False


def resolve_records(batch, callermap):
    if len(batch) == 1:
        for r in resolve_merge(batch, callermap):
            yield r
    else:
        logger.info("Processing batch size: %s", len(batch))
        if contains_pindel(batch):
            pc = [r for r in batch if 'pindel' in r.caller][0]
            callset = []
            for r in batch:
                if 'pindel' in r.caller:
                    callset.append(callermap[r.caller])
                else:
                    callset.append(callermap[r.caller] + '*')
            logger.info("Merged pindel call with %s", callset)
            yield resolve_merge([pc], callermap)[0], '|'.join(callset)
        else:
            grouped = defaultdict(list)
            for r in batch:
                for alt in r.record['ALT']:
                    key = (r.record['CHROM'], r.record['POS'], r.record['REF'], alt)
                    new_record = r.record.copy()
                    new_record['ALT'] = [alt]
                    grouped[key].append(MetaRecord(r.caller, new_record))

            for key, records in grouped.items():
                for r in resolve_merge(records, callermap):
                    yield r


def main(args):
    reader = MultiVCFReader(args.input, args.output, args.keys)
    reader.outwriter.header.add_header('##INFO=<ID=CENTERS,Number=1,Type=String,Description="Center files that made the call">')
    reader.outwriter.write_header()
    
    callermap = dict(zip(args.input, args.keys))

    for chunk in reader:
        for r, c in resolve_records(chunk, callermap):
            r['INFO'] = {'CENTERS':[c]}
            reader.outwriter.write_record(r)
    logger.info("Done")
            

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument('--keys', type = str, help = 'caller keys', nargs = '+')
    parser.add_argument('--output', type = str, help = 'output file')
    parser.add_argument('--input', nargs='+', type = str, help = 'input files')

    args = parser.parse_args()
    main(args)

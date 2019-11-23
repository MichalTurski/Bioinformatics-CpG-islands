import pandas as pd
import spans

#%% read files
islands_df = pd.read_csv('cpgIslandExt.txt', sep='\t', usecols=[1, 2, 3], header=None)
islands_df.columns=['chr', 'start', 'stop']
lengths_df = pd.read_csv('hg19.chrom.sizes.txt', sep='\t', index_col=0, header=None, squeeze=True)
lengths = lengths_df.to_dict()
methyalation_df = pd.read_csv('HAIB.A549.EtOH.Rep.3.bed', sep='\t', usecols=[0, 1, 2], header=None)
methyalation_df.columns=['chr', 'start', 'stop']

#%% filter inputs
autosomal_chrom = set()
for i in range(22):
    autosomal_chrom.add('chr' + str(i + 1))


def filter_autosomal(df, autosomal_chrom):
    return df[df['chr'].isin(autosomal_chrom)]


islands_df = filter_autosomal(islands_df, autosomal_chrom)
# lengths_df = filter_autosomal(lengths_df, autosomal_chrom)
methyalation_df = filter_autosomal(methyalation_df, autosomal_chrom)

#%% find islands, shores, shelves and seas


def intervals_for_chrom(islands_chrom_df, chrom_len):
    chrom_interval = spans.intrange(0, chrom_len)
    sea_interval = spans.intrangeset(chrom_interval)
    shelf_interval = spans.intrangeset(None)
    shore_interval = spans.intrangeset(None)
    island_interval = spans.intrangeset(None)
    for _, island in islands_chrom_df.iterrows():
        island_interval.add(spans.intrange(island['start'], island['stop']).intersection(chrom_interval))


for chrom in autosomal_chrom:
    spans.intrange(0, lengths[chrom])
    sea_interval = spans.intrangeset([spans.intrange(0, lengths[chrom])])
    shelf_interval = spans.intrangeset([])
    shore_interval = spans.intrangeset([])



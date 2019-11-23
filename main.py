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
for i in range(3):  # TODO: range to 22
    autosomal_chrom.add('chr' + str(i + 1))


def filter_autosomal(df, autosomal_chrom):
    return df[df['chr'].isin(autosomal_chrom)]


islands_df = filter_autosomal(islands_df, autosomal_chrom)
# lengths_df = filter_autosomal(lengths_df, autosomal_chrom)
methyalation_df = filter_autosomal(methyalation_df, autosomal_chrom)

#%% find islands, shores, shelves and seas


def intervals_for_chrom(islands_chrom_df, chrom_len):
    # Border numbers section:
    shore_limit = 2000
    shelf_limit =  shore_limit + 2000

    chrom_interval = spans.intrange(0, chrom_len)
    intervals = {'sea': spans.intrangeset([chrom_interval]), 'shelf': spans.intrangeset([]),
                 'shore': spans.intrangeset([]), 'island': spans.intrangeset([])}
    for _, island in islands_chrom_df.iterrows():
        intervals['island'].add(spans.intrange(island['start'], island['stop']).intersection(chrom_interval))
        intervals['shore'].add(spans.intrange(island['start'] - shore_limit, island['stop'] + shore_limit)
                               .intersection(chrom_interval))
        intervals['shelf'].add(spans.intrange(island['start'] - shelf_limit, island['stop'] + shelf_limit)
                               .intersection(chrom_interval))
    print('done inner loop')
    intervals['sea'] = intervals['sea'].difference(intervals['shelf'])
    print('done sea interval')
    intervals['shelf'] = intervals['shelf'].difference(intervals['shore'])
    print('done shelf interval')
    intervals['shore'] = intervals['shore'].difference(intervals['island'])
    print('done shore interval')

    return intervals


def df_from_interval(interval, chrom):
    fragments_list =[]
    for subinterval in interval:
        fragments_list.append((chrom, subinterval.lower, subinterval.upper))
    df = pd.DataFrame(fragments_list)
    return df


dfs = {'island': [], 'shore': [], 'shelf': [], 'sea': []}
chrom_interval_dict = {}
for chrom in autosomal_chrom:
    islands_chrom_df = islands_df[islands_df['chr'] == chrom]
    intervals = intervals_for_chrom(islands_chrom_df, lengths[chrom])
    for key in intervals:
        dfs[key].append(df_from_interval(intervals[key], chrom))
    chrom_interval_dict[chrom] = intervals
    # dfs['island'].append(df_from_interval(island_interval, chrom))
    # dfs['shore'].append(df_from_interval(shore_interval, chrom))
    # dfs['shelf'].append(df_from_interval(shelf_interval, chrom))
    # dfs['sea'].append(df_from_interval(sea_interval, chrom))


new_islands_df = pd.concat(dfs['island'])
shores_df = pd.concat(dfs['shore'])
shelves_df = pd.concat(dfs['shelf'])
seas_df = pd.concat(dfs['sea'])

#%% save results to files
new_islands_df.to_csv('islands.bed', sep='\t', header=False, index=False)
shores_df.to_csv('shores.bed', sep='\t', header=False, index=False)
shelves_df.to_csv('shelves.bed', sep='\t', header=False, index=False)
seas_df.to_csv('seas.bed', sep='\t', header=False, index=False)

#%% calculate methylations positions on chromosomes
methyalation_df['pos'] = (methyalation_df['stop'] - methyalation_df['start'])/2 + methyalation_df['start']
methyalation_df['pos'] = methyalation_df['pos'].apply(lambda x: int(x))

#%% calculate methylations locations in islands/shores/shelves/seas
def detemine_methylation_location(chrom, pos, chrom_interval_dict):
    intervals = chrom_interval_dict[chrom]
    for region in intervals:
        if intervals[region].contains(pos):
            return region
    raise LookupError(f'There is no position {str(pos)} in chromosome {chrom}.')


methyalation_df['location'] = methyalation_df.apply(
    lambda row: detemine_methylation_location(row['chr'], row['pos'], chrom_interval_dict), axis=1)

#%%





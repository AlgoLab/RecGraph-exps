import numpy as np

def table(times_ab, rams_ab, times_rs, rams_rs , out='md'):
    if out == 'md':
        print('| | Average Time (s) | std Time (s) | Average RAM (GB) | std RAM (GB) |')
        print('| --- | ---: | ---: | ---: | ---: |')
        print(f'| abPOA | {np.mean(times_ab):.4f} | {np.std(times_ab):.4f} | {np.mean(rams_ab):.4f} | {np.std(rams_ab):.4f} |')
        print(f'| rsPOA | {np.mean(times_rs):.4f} | {np.std(times_rs):.4f} | {np.mean(rams_rs):.4f} | {np.std(rams_rs):.4f} |')

def describe(_times, _rams):
    print(f'count\t{len(_times)}')
    print(f'mean\t\t{np.mean(_times):10.4f} s\t{np.std(_rams):10.5f} GB')
    print(f'std\t\t{np.std(_times):10.4f} s\t{np.std(_rams):10.5f} GB')
    print(f'min\t\t{min(_times):10.4f} s\t{min(_rams):10.5f} GB')
    print(f'max\t\t{max(_times):10.4f} s\t{max(_rams):10.5f} GB')
    print(f'sum\t\t{sum(_times):10.4f} s\t{sum(_rams):10.5f} GB')



def load_data(files, tool):
    times = []
    rams = []
    avg_trun = []
    avg_rrun = []

    for f in files:
        tt = []
        rr = []
        for line in open(f, 'r'):
            if tool == 'abpoa':
                if '[abpoa_main]' in line:
                    data = line.strip().split(' ')
                    tt.append(float(data[3]))
                    rr.append(float(data[10]))
            if tool == 'rspoa':
                if 'Elapsed (wall clock)' in line:
                    t = line.strip().split(' ')[-1]
                    m, sec = t.split(':')
                    tt.append(int(m)*60 + float(sec))
                elif 'Maximum resident set' in line:
                    r = line.strip().split(' ')[-1]
                    rr.append(float(r)*1e-6) # TODO: CHECKME

        avg_trun.append(np.mean(tt))
        times += tt
        avg_rrun.append(np.mean(rr))
        rams += rr

    
    return times, rams, avg_trun, avg_rrun

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('--abpoa', type=str, nargs='+',
                        help='abpoa times file')
    parser.add_argument('--rspoa', type=str, nargs='+',
                        help='rspoa times file')

    args = parser.parse_args()

    ab_times, ab_rams, ab_avg_trun, ab_avg_rrun = load_data(args.abpoa, tool='abpoa')
    print( '-- abPOA ' + '-'*30)
    describe(ab_times, ab_rams)
    
    rs_times, rs_rams, rs_avg_trun, rs_avg_rrun = load_data(args.rspoa, tool='rspoa')
    print( '-- rsPOA ' + '-'*30)
    describe(rs_times, rs_rams)
    

    print( '-- MD Table ' + '-'*30)
    table(ab_times, ab_rams, rs_times, rs_rams)
    print()
    print(f'**abPOA** avg per region: \t{np.mean(ab_avg_trun):.4f} s, \t{np.mean(ab_avg_rrun):.4f} GB \\')
    print(f'**rsPOA** avg per region: \t{np.mean(rs_avg_trun):.4f} s, \t{np.mean(rs_avg_rrun):.4f} GB')
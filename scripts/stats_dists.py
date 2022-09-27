import numpy as np

def table(dists_ab, dists_rs, out='md'):
    if out == 'md':
        print('| | Average Dist | std Dist | Max Dist | Min Dist |')
        print('| --- | ---: | ---: | ---: | ---: |')
        if len(dists_ab) > 0:
            print(f'| abPOA | {np.mean(dists_ab):.2f} | {np.std(dists_ab):.2f} | {max(dists_ab):.2f} | {min(dists_ab):.2f} |')
        if len(dists_rs) > 0:
            print(f'| rsPOA | {np.mean(dists_rs):.2f} | {np.std(dists_rs):.2f} | {max(dists_rs):.2f} | {min(dists_rs):.2f} |')

def describe(_dists):
    print(f'count\t{len(_dists)}')
    print(f'mean\t\t{np.mean(_dists):10.4f}')
    print(f'std\t\t{np.std(_dists):10.4f}')
    print(f'min\t\t{min(_dists):10.4f}')
    print(f'max\t\t{max(_dists):10.4f}')

def load_data_percalns(files):
    percs = []
    avg_prun = []

    if not files:
        return percs, avg_prun

    for f in files:
        pp = []
        for line in open(f, 'r'):
            pp.append(float(line.strip()))

        avg_prun.append(np.mean(pp))
        percs += pp

    return percs, avg_prun

def load_data_dists(files):
    dists = []
    avg_drun = []

    if not files:
        return dists, avg_drun

    for f in files:
        dd = []
        for line in open(f, 'r'):
            dd.append(float(line.strip()))

        avg_drun.append(np.mean(dd))
        dists += dd

    return dists, avg_drun

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('--abpoa', type=str, nargs='+',
                        help='abpoa levdist file(s)')
    parser.add_argument('--rspoa', type=str, nargs='+',
                        help='rspoa levdist file(s)')
    parser.add_argument('--rspercal', type=str, nargs='+',
                        help='rspoa percalns file(s)')
    

    args = parser.parse_args()
  
    ab_dists, ab_avg_drun = load_data_dists(args.abpoa)
    if args.abpoa:
        print( '-- abPOA (edit dist)' + '-'*30)
        describe(ab_dists)

    rs_dists, rs_avg_drun = load_data_dists(args.rspoa)
    if args.rspoa:
        print( '-- rsPOA (edit dist)' + '-'*30)
        describe(rs_dists)

    table(ab_dists, rs_dists)
    print()
    if args.abpoa:
        print(f'**abPOA** avg per region: \t{np.mean(ab_avg_drun):.4f} \\')
    if args.rspoa:
        print(f'**rsPOA** avg per region: \t{np.mean(rs_avg_drun):.4f} \\')

    rs_pa, rs_avg_prun = load_data_dists(args.rspercal)
    if args.rspercal:
        print(f'**rsPOA** avg % of read alignment: \t{np.mean(rs_pa):.4f}')
        # print( '-- rsPOA (%% read alignment)' + '-'*30)
        # describe(rs_pa)
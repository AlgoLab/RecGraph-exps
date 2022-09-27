import numpy as np

def table(dists_ab, dists_rs, out='md'):
    if out == 'md':
        print('| | Average Dist | std Dist | Max Dist | Min Dist |')
        print('| --- | ---: | ---: | ---: | ---: |')
        print(f'| abPOA | {np.mean(dists_ab):.2f} | {np.std(dists_ab):.2f} | {max(dists_ab):.2f} | {min(dists_ab):.2f} |')
        print(f'| rsPOA | {np.mean(dists_rs):.2f} | {np.std(dists_rs):.2f} | {max(dists_rs):.2f} | {min(dists_rs):.2f} |')

def describe(_dists):
    print(f'count\t{len(_dists)}')
    print(f'mean\t\t{np.mean(_dists):10.4f}')
    print(f'std\t\t{np.std(_dists):10.4f}')
    print(f'min\t\t{min(_dists):10.4f}')
    print(f'max\t\t{max(_dists):10.4f}')

def load_data(files):
    dists = []
    avg_drun = []

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

    args = parser.parse_args()
    
    ab_dists, ab_avg_drun = load_data(args.abpoa)
    print( '-- abPOA ' + '-'*30)
    describe(ab_dists)

    rs_dists, rs_avg_drun = load_data(args.rspoa)
    print( '-- rsPOA ' + '-'*30)
    describe(rs_dists)

    table(ab_dists, rs_dists)
    print()
    print(f'**abPOA** avg per region: \t{np.mean(ab_avg_drun):.4f} \\')
    print(f'**rsPOA** avg per region: \t{np.mean(rs_avg_drun):.4f}')
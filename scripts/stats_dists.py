import numpy as np

def describe(_dists):
    print(f'count\t{len(_dists)}')
    print(f'mean\t\t{np.mean(_dists):10.4f}')
    print(f'std\t\t{np.std(_dists):10.4f}')
    print(f'min\t\t{min(_dists):10.4f}')
    print(f'max\t\t{max(_dists):10.4f}')
    print(f'sum\t\t{sum(_dists):10.4f}')



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
    parser.add_argument('--rspoa', type=str, nargs='+',
                        help='rspoa levdist file(s)')

    args = parser.parse_args()
    
    rs_dists, rs_avg_drun = load_data(args.rspoa)
    print( '-- rsPOA ' + '-'*30)
    describe(rs_dists)
    print()
    print(f'**rsPOA** avg per region: \t{np.mean(rs_avg_drun):.4f}')
import requests

def download_abtem(directory, verbose=False, debug=False):
    """
    This function downloads parametrization data files from abTEM.
    """
    files = ('kirkland.json',
             'lobato.json',
             'peng_high.json',
             'peng_ionic.json',
             'peng_low.json',
             'waasmaier_kirfel.json')

    with requests.Session() as s:
        for file in files:
            url = 'https://raw.githubusercontent.com/abTEM/abTEM/' + \
                f'main/abtem/parametrizations/data/{file}'
                
            if verbose or debug:
                print(f'Downloading {url} ...')
            r = s.get(url)
            
            if verbose or debug:
                print('Success!\n')
            
            if debug:
                print(r.headers.keys())
                print(r.__dict__.keys())
                print(r._content)
            
            
            ofile = f'{directory}/{file}'
            if verbose or debug:
                print(f'Writing data file {ofile} ...')
            
            with open(ofile,'w') as f:
                f.write(r.text)
                f.write('\n')
            
            if verbose or debug:
                print('Success!\n')

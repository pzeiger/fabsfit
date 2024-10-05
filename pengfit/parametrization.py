import json
from scipy.optimize import curve_fit

class Parametrization():

    def __init__(self, ):
        pass
    
    def fit(self):
        popt, pcov = curve_fit(self.scattering_factor, xdata, ydata)
        pass
    
    def load_data(self):
        self.data = json.load(self._datafile)
    
    def scattering_factor(self):
        pass


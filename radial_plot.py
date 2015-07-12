def main():

    import h5py
    import numpy
    import matplotlib.pyplot as plt

    with h5py.File('final.h5','r+') as f:
        raw = {}
        for field in ['x_coordinate',
                      'y_coordinate',
                      'density',
                      'pressure',
                      'x_velocity',
                      'y_velocity',
                      'ghost']:
            raw[field] = numpy.array(f[field])
        raw['radius'] = numpy.sqrt(raw['x_coordinate']**2+
                                   raw['y_coordinate']**2)
        raw['r_velocity'] = ((
            raw['x_velocity']*raw['x_coordinate']+
            raw['y_velocity']*raw['y_coordinate'])/
                             raw['radius'])
        mask = (raw['ghost']<0.5)
        for n,field in enumerate(['density','pressure','r_velocity']):
            plt.subplot(3,1,n+1)
            plt.plot(raw['radius'][mask],
                     raw[field][mask],
                     '.')
        plt.show()

if __name__ == '__main__':

    main()

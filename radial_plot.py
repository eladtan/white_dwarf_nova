def main():

    import h5py
    import numpy
    import matplotlib.pyplot as plt

    fi = h5py.File('initial.h5','r+')
    with h5py.File('final.h5','r+') as f:
        raw = {}
        for field in ['x_coordinate',
                      'y_coordinate',
                      'density',
                      'temperature',
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
        for n,field in enumerate(['density','pressure','temperature','y_velocity']):
            plt.subplot(2,2,n+1)
            initial_radius = numpy.sqrt(
                numpy.array(fi['x_coordinate'])**2+
                numpy.array(fi['y_coordinate'])**2)
            plt.plot(initial_radius[mask],
                     numpy.array(fi[field])[mask],
                     '.')
            plt.plot(raw['radius'][mask],
                     raw[field][mask],
                     '.')
            plt.ylabel(field)
        plt.show()

if __name__ == '__main__':

    main()

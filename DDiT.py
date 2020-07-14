import sys
import time
import numpy as np
import matplotlib.pyplot as plt
# -----------------------------------------------------------
np.seterr(over='ignore')
# -----------------------------------------------------------
class Disk(object):
    """
    docstring for RT
    """
    def __init__(self, nx = 300, pixscale = 0.01226, nframe = 50):
        """
        description
        """
        if type(nx) is not int: self._error_msg('The parameter \'nx\' should be an integer ')
        if type(nframe) is not int: self._error_msg('The parameter \'nframe\' should be an integer ')
        if nframe%2 == 1:
            print('It is probably better if \'nframe\' is an even number, to sample to midplace properly.')
        self._threshold = 5.e-2
        self._nx = nx
        self._cx = self._nx//2
        self._nframe = nframe
        self._pixscale = pixscale
        self._xlim = self._cx * self._pixscale
        self._Xin, self._Yin = (np.mgrid[0:self._nx, 0:self._nx] - self._cx) * self._pixscale
        self._xm, self._ym = np.zeros(shape=(self._nx, self._nx)), np.zeros(shape=(self._nx, self._nx))
        """
        Define some variables that will contain the:
            - total intensity image
            - polarized intensity image
            - scattering angle at the midplane
            - distance to the star at the midplane
            - azimuth angle at the midplane
        """
        self.polarized = np.zeros(shape=(self._nx, self._nx))
        self.intensity = np.zeros(shape=(self._nx, self._nx))
        self.scattering = np.zeros(shape=(self._nx, self._nx))
        self.distance = np.zeros(shape=(self._nx, self._nx))
        self.azimuth = np.zeros(shape=(self._nx, self._nx))
        """
        Some trigonometry variables
        """
        self._cs, self._ss, self._st, self._st2, self._to = 0., 0., 0., 0., 0.
        """
        The geometric parameters for the disk
        """
        self.a, self.incl, self.pa, self.pin, self.pout, self.e, self.omega, self.opang = 1., 0.1, 152.1*np.pi/180., 25., -2.5, 0., 0.,0.04
        """
        Parameters for the phase functions
        """
        self.gsca, self.gpol = 0., 0.
        self.theta, self.s11, self.s12 = None, None, None

        self._ismodel = False

    """
    Get the entry and exit points of an ellipse in the [y,z] plane
    """
    def _get_sphere(self):
        """
        Intersection points with an ellipse along the [y,z] plane
        """
        ye, ys = np.zeros(shape=(self._nx, self._nx)), np.zeros(shape=(self._nx, self._nx))
        ze, zs = np.zeros(shape=(self._nx, self._nx)), np.zeros(shape=(self._nx, self._nx))
        xe = np.zeros(shape=(self._nx, self._nx))
        """
        it should be a * (1 - e**2) / (1-e) which simplyfies as a*(1+e)
        """
        rmax = self._threshold**(1.e0 / self.pout) * self.a * (1.e0 + self.e) 
        #print(rmax)
        radius = rmax**2. - self._xm**2.
        radius[(radius<0)] = 0.
        radius = np.sqrt(radius)
        height = 3. * radius * self._to

        delta = self._ym**2. * self._st2**2. * radius ** 4. - (height**2. + self._st2 * radius**2.) * (self._st2 * self._ym**2. * radius**2. - height**2. * radius**2.)
        sel = (delta > 0.) 
        ye[sel] = ((self._ym[sel] * self._st2 * radius[sel]**2.) + np.sqrt(delta[sel])) / (height[sel]**2. + self._st2 * radius[sel]**2.)
        ys[sel] = ((self._ym[sel] * self._st2 * radius[sel]**2.) - np.sqrt(delta[sel])) / (height[sel]**2. + self._st2 * radius[sel]**2.)
        ze[sel] = self._st * (ye[sel] - self._ym[sel])
        zs[sel] = self._st * (ys[sel] - self._ym[sel])
        xe[sel] = self._xm[sel]

        return xe, ye, ze, ys, zs 

    """
    Compute the model for the total and polarized intensity.
    """
    def compute_model(self, **kwargs):
        """
        Method to compute a synthetic image of a debris disk

        INPUTS:
          - a: reference semi-major axis in arcseconds
          - incl: inclination in degrees
          - pa: position angle in degrees
          - pin: inner slope
          - pout: outer slope
          - gsca: HG coefficient for total intensity
          - gpola: HG coefficient for polarized intensity
          - e: eccentricity
          - omega: argument of pericenter in degrees
          - opang: opening angle of the disk
          - s11: the scattered light phase function 
          - s12: the polarized light phase function 
        """
        """
        Check the input parameters
        And put the angles in radians
        """
        self._check_parameters(kwargs)
        """
        Compute some trigonometric things
        """
        self._trigonometry()
        """
        Need to re-initialize the arrays to zero, in case I compute two models one after the other
        """
        self.polarized = np.zeros(shape=(self._nx, self._nx))
        self.intensity = np.zeros(shape=(self._nx, self._nx))
        """
        Project the disk according to the inclination and position angle
        """
        self._xm = ((np.cos(self.pa) * self._Xin + np.sin(self.pa) * self._Yin))
        self._ym = ((np.sin(self.pa) * self._Xin - np.cos(self.pa) * self._Yin)) / self._cs
        self.distance = np.sqrt(self._xm**2. + self._ym**2.)
        self.distance[self._cx, self._cx] = 1.
        self.scattering = np.arccos((self._ss * self._ym)/self.distance)
        """
        Get the entry and exit points for the upper and lower layers
        """
        xe, ye, ze, ys, zs = self._get_sphere()
        self._get_flux(xe, ye, ze, ys, zs)
        self._ismodel = True

    """
    Method to compute the emission between the entry and exit points.
    """
    def _get_flux(self, xe, ye, ze, ys, zs):
        """
        Compute the azimuth angle in the disk midplane
        This will put the north side at 
        """
        self.azimuth = (np.arctan2(self._ym,self._xm) + np.pi/2.) % (2. * np.pi) - np.pi
        """
        The volume of each cells
        The line of sight is divided as follows:
        for i in range(nframe):
            zi = (ze-zs) * i / (nframe - 1) + zs
        So then I can compute the lenght of a given cell which is (ze-zs)/(nframe -1)
        And the volume of one cell is sqrt(Delta_z**2. + Delta_y**2) * pixelscale**2.
        """
        volume = np.sqrt(((ze-zs)**2. + (ye-ys)**2.))/(self._nframe-1.)*self._pixscale**2.
        for i in range(self._nframe-1):
            """
            Get the middle points of each cells.
            """
            zi = zs + (ze - zs) * (i+0.5) / (self._nframe-1)
            yi = ys + (ye - ys) * (i+0.5) / (self._nframe-1)
            """
            Define variables for the phase function, that need to be re-initialized for every iteration
            This slows down the code a bit, but it is more accurate and prevents some unfortunate 
            NaNs in some cases (especially if the position angle is 0 or 90 degrees.

            I had tried to define them once, and set them to 0., but that did not improve the time.
            For clarity, I leave it as it is now.
            """
            psca = np.zeros(shape=(self._nx, self._nx, 2))
            ppol = np.zeros(shape=(self._nx, self._nx, 2))
            theta = np.zeros(shape=(self._nx, self._nx))
            costheta = np.zeros(shape=(self._nx, self._nx))
            densr = np.zeros(shape=(self._nx, self._nx))
            densz = np.zeros(shape=(self._nx, self._nx))
            density = np.zeros(shape=(self._nx, self._nx))
            """
            To avoid some issues for very inclined disks, I need to re-evaluate
            the azimuthal angles for each "frames", as well as the reference
            radius, based on the azimuthal angle.

            This azimuthal angle is the one in the midlpane because I need to compute the 
            reference radius at the midplane.
            """
            az = np.arctan2(yi, xe)
            r_ref = self.a * (1.e0 - self.e * self.e) / (1.e0 + self.e * np.cos(az + self.omega))
            """
            Compute the cosine of the scattering angle
            I need to make a selection where dist3d != 0, to avoid division by 0 when computing it.
            """
            dist3d = np.sqrt(xe**2. + yi**2. + zi**2.)
            sel3d = (dist3d > 0.)
            if self.theta is None:
                costheta[sel3d] = (self._ss * yi[sel3d] - self._cs * zi[sel3d]) / dist3d[sel3d]
                psca[sel3d,0] = (1.e0 - self.gsca**2.) / (4. * np.pi * (1.e0 + self.gsca**2. - 2. * self.gsca * costheta[sel3d])**(1.5))
                psca[sel3d,1] = (1.e0 - self.gsca**2.) / (4. * np.pi * (1.e0 + self.gsca**2. - 2. * self.gsca * costheta[sel3d])**(1.5))
                ppol[sel3d,0] = (1.e0 - self.gpol**2.) / (4. * np.pi * (1.e0 + self.gpol**2. - 2. * self.gpol * costheta[sel3d])**(1.5)) * (1.e0 - costheta[sel3d]**2.) / (1.e0 + costheta[sel3d]**2.)
                ppol[sel3d,1] = (1.e0 - self.gpol**2.) / (4. * np.pi * (1.e0 + self.gpol**2. - 2. * self.gpol * costheta[sel3d])**(1.5)) * (1.e0 - costheta[sel3d]**2.) / (1.e0 + costheta[sel3d]**2.)
            else:
                theta[sel3d] = np.arccos((self._ss * yi[sel3d] - self._cs * zi[sel3d]) / dist3d[sel3d])
                psca[sel3d,0] = np.interp(theta[sel3d], self.theta, self.s11[:,0])
                psca[sel3d,1] = np.interp(theta[sel3d], self.theta, self.s11[:,1])
                ppol[sel3d,0] = np.interp(theta[sel3d], self.theta, self.s12[:,0])
                ppol[sel3d,1] = np.interp(theta[sel3d], self.theta, self.s12[:,1])
            """
            Define the density structure
            I also need to check if dist2d, the distance in the midplane, is null or not. If it is 0
            then it may yield to some problems in the exponent when computing the radial density profile.
            """
            dist2d = np.sqrt(xe**2. + yi**2.)
            sel2d = (dist2d > 0.)
            densr[sel2d] = ((dist2d[sel2d]/r_ref[sel2d])**(-2.*self.pout) + (dist2d[sel2d]/r_ref[sel2d])**(-2.*self.pin))**(-.5)
            densz[sel2d] = np.exp(-zi[sel2d]**2 / (2. * (self._to * dist2d[sel2d])**2.))
            density[sel3d] = densr[sel3d] * densz[sel3d] * volume[sel3d] / (dist3d[sel3d]**2.)
            """
            Differentiate the north and south sides of the disk
            North side is where azi_tmp is negative.

            When doing:
            azi_tmp = (self.azimuth + np.pi/2.) % (2. * np.pi) - np.pi
            This should put the angle at 0. along the minor axis of the disk and between -pi and +pi for both sides of the disk.
            """
            sel = (self.azimuth <=0.) 
            self.intensity[sel] += density[sel] * psca[sel,0]
            self.intensity[~sel] += density[~sel] * psca[~sel,1]
            self.polarized[sel] += density[sel] * ppol[sel,0]
            self.polarized[~sel] += density[~sel] * ppol[~sel,1]

            #tmp = np.zeros(shape=(self._nx, self._nx))
            #tmp[sel] = density[sel] * psca[sel,0]
            #tmp[~sel] = density[~sel] * psca[~sel,1]

            #xlim = self._cx * self._pixscale
            #fig = plt.figure(figsize=(7,7 * 9./16.))
            #ax1 = fig.add_axes([0.0, 0.0, 1., 1.])
            ##ax1.imshow(tmp, origin = 'lower', vmin = 0., vmax = np.pi, cmap = 'inferno')
            #ax1.imshow(tmp, origin = 'lower', vmin = 0., vmax = 3.5e-7, cmap = 'inferno', extent=[xlim, -xlim, -xlim, xlim])
            #ax1.plot(0.,0., marker = '+', ms = 6 , color = 'w')
            #ax1.set_xlim(xlim, -xlim)
            #ax1.set_ylim(-xlim * 9./16., xlim * 9./16.)
            #ax1.axis('off')
            #plt.savefig('debug/image_'+format(i,'04d')+'.png', edgecolor='black', dpi=100, facecolor='black')
            ##if i==0:
                ##plt.show()
            #plt.close()

    """
    Plot the images
    """
    def plot(self, cmap = 'inferno'):
        """
        Method to make some plot
        """
        if self._ismodel:
            fig = plt.figure(figsize=(13,6))
            ax1 = fig.add_axes([0.1, 0.12, 0.4, 0.8])
            im = ax1.imshow(self.intensity, origin = 'lower', extent = [self._xlim, -self._xlim, -self._xlim, self._xlim], cmap = cmap, vmin = np.percentile(self.intensity, 1.), vmax = np.percentile(test.intensity, 99.9))
            ax1.set_xlabel('$\Delta \\alpha$ [$^{\prime\prime}$]')
            ax1.set_ylabel('$\Delta \delta$ [$^{\prime\prime}$]')

            ax2 = fig.add_axes([0.55, 0.12, 0.4, 0.8])
            im = ax2.imshow(self.polarized, origin = 'lower', extent = [self._xlim, -self._xlim, -self._xlim, self._xlim], cmap = cmap, vmin = np.percentile(self.polarized, 1.), vmax = np.percentile(test.polarized, 99.9))
            ax2.set_xlabel('$\Delta \\alpha$ [$^{\prime\prime}$]')
            plt.show()

    """
    Check the parameters that are passed as kwargs
    """
    def _check_parameters(self, kwargs):
        """
        Check parameters that are being passed and make 
        some consistency checks for the inclination
        """
        if 'a' in kwargs:
            self.a = kwargs['a']
        if 'incl' in kwargs:
            self.incl = kwargs['incl'] * np.pi / 180.
        if 'pa' in kwargs:
            self.pa = -kwargs['pa'] * np.pi / 180.
        if 'pin' in kwargs:
            self.pin = kwargs['pin']
        if 'pout' in kwargs:
            self.pout = kwargs['pout']
        if 'gsca' in kwargs:
            self.gsca = kwargs['gsca']
        if 'gpol' in kwargs:
            self.gpol = kwargs['gpol']
        if 'e' in kwargs:
            self.e = kwargs['e']
        if 'omega' in kwargs:
            self.omega = kwargs['omega'] * np.pi / 180.
        if 'opang' in kwargs:
            self.opang = kwargs['opang']
        if 's11' in kwargs:
            self.s11 = kwargs['s11']
            self.theta = np.linspace(0., np.pi, num = np.shape(self.s11)[0])
            if 's12' not in kwargs: self.s12 = np.ones(shape=(np.shape(self.s11)[0],2))
        if 's12' in kwargs:
            self.s12 = kwargs['s12']
            self.theta = np.linspace(0., np.pi, num = np.shape(self.s12)[0])
            if 's11' not in kwargs: self.s11 = np.ones(shape=(np.shape(self.s12)[0],2))

    """
    Define some cos, sin, and tan
    """
    def _trigonometry(self):
        """
        Define some variables
        """
        if self.incl == 0.:
            self.incl = 0.1 * np.pi / 180.
        if self.incl == np.pi/2.:
            self.incl = 89.9 * np.pi / 180.
        self._cs, self._ss = np.cos(self.incl), np.sin(self.incl)
        self._ti = np.tan(self.incl)
        self._st = self._cs / self._ss
        self._st2 = self._st**2.
        self._to = np.tan(self.opang)

    """
    Print some error message and quit
    """
    def _error_msg(self, message):
        """
        Print something and quit
        """
        print(message)
        sys.exit()

if __name__ == '__main__':        
    test = Disk(nframe = 50)
    t0 = time.time()
    test.compute_model(e = 0.1, incl = 69., pa = 110., a = 0.89, gsca = 0.4, gpol = 0.6, omega = 180., opang = 0.05, pin = 20.0, pout = -5.5)
    print('Took: ' + format(time.time()-t0, '0.2f') + ' seconds.')
    test.plot()


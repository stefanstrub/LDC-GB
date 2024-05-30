""" Project waveforms h+ and hx on LISA arms. """

import h5py
import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline as spline
import lisaconstants
import ldc.io.hdf5 as hdfio
try:
    import cupy as cp
except ImportError:
    cp = None

#pylint:disable=C0103
#pylint:disable=W0201

clight = lisaconstants.SPEED_OF_LIGHT
half_schwarzschild_radius = lisaconstants.SUN_SCHWARZSCHILD_RADIUS/2
post_newtonian_constant = 1.

LIST_SEP = " , "

def from_file(hdf5_filename, nodata=False, ilink=None):
    """Load projected strain from hdf5 file

    if nodata is True, then only meta data are read.  ilink can be
    used to load a single link (index between 0 and 5)
    """
    if nodata:
        attrs = hdfio.load_attributes(hdf5_filename, name="strain")
        yArm = None
    else:
        yArm, attrs = hdfio.load_array(hdf5_filename, name="strain")
    source_names = attrs['source_names'].split(LIST_SEP)
    links = attrs['links']
    links = links.split(LIST_SEP)
    t_min =  attrs['t_min']
    t_max = attrs["t_max"]
    dt = attrs["dt"]
    if ilink is not None:
        yArm = yArm[:,ilink]
    return yArm, source_names, links, t_min, t_max, dt

def to_file(hdf5_filename, yArm, source_names, links, t_min, t_max, dt, ilink=None):
    """Save projected strain into hdf5 file.

    Meta data are:

        - list of source names summed in the strain
        - list of link name like (1-2), (2-3) etc. for (receiver-emitter)
        - time vector description : t_min, t_max, dt

    ilink can be used to write a single link at a time. Meta data are
    written for ilink=0.

    """
    if ilink is None or ilink==0:
        if len(source_names)>100:
            str_name = "-".join([source_names[0], source_names[-1]])
        else:
            str_name = LIST_SEP.join(source_names)
        str_link = LIST_SEP.join(links)
        chunks = True if ilink==0 else False
        hdfio.save_array(hdf5_filename, yArm, name='strain', mode='w',
                         chunks=chunks, links=str_link,
                         source_names=str_name, t_min=t_min, t_max=t_max, dt=dt)
    else:
        hdfio.append_array(hdf5_filename, yArm, ilink, name='strain')


class ProjectedStrain():
    """ Project GW strain on LISA arms """

    def __init__(self, orbits):
        """  Set orbits and constellation configuration.
        """
        self.orbits = orbits
        self.set_link()

    def init_links(self, receiver_time, order=0, cuda=False):
        """ Compute and save sc positions and travel time.
        """

        if cuda is True and cp is None:
            print("cupy is not installed, can't use cuda option")
            cuda = False
        #alphas = self.orbits.compute_alpha(receiver_time)

        # Choose between numpy and cupy
        xp = cp if cuda is True else np

        # pos and tt will be saved as numpy/cupy arrays
        self.pos = xp.zeros((self.orbits.number_of_spacecraft, 3, len(receiver_time)))
        for i in range(1,self.orbits.number_of_spacecraft+1):
            self.pos[i-1,:,:] = xp.array(self.orbits.compute_position(i,
                                                                      receiver_time)/clight)

        self.tt = np.zeros((len(receiver_time), self.orbits.number_of_arms))
        for i in range(self.orbits.number_of_arms):
            receiver, emitter = self.orbits.get_pairs()[i]
            pe = self.pos[emitter-1,:]*clight
            pr = self.pos[receiver-1,:]*clight
            self.tt[:,i] = self.orbits.compute_travel_time(emitter, receiver,
                                                           receiver_time, order=order,
                                                           position_emitter=pe,
                                                           position_receiver=pr)


    def arm_response(self, t_min, t_max, dt, GWs, tt_order=0, #extended_t=1e3,
                     interp_type=['EMRI', 'MBHB', 'Numeric'], ilink=None, **kwargs):
        """Return projected strain y_t

        For source in interp_type, hphc is computed only once, then
        stored and used for interpolation. Beware of the memory
        consumption of such sources.

        Args:
            t_min, t_max, dt: scalar defining time range
            GWs: list of HpHc object
            ilink: specify the link number (0 to 5) on which to project the signal.
                   If None, all links are computed.

        >>> proj = ProjectedStrain(orbits)
        >>> yArm = proj.arm_response(0, 1e4, 10, [GB])
        >>> print(yArm[0,:])
        [-5.83708114e-24  4.67226745e-24 -5.83693542e-24 -1.03547699e-25
          4.67357629e-24 -1.03579292e-25]
        """
        receiver_time = np.arange(t_min, t_max, dt)
        #extended_time = np.arange(max(t_min-dt*(extended_t//dt), 0),
        #                              t_max+dt*(extended_t//dt), dt)

        ## GW hp,hc precomputation
        self.source_names = [GW.source_name for GW in GWs]
        #jk = [GW.compute_hphc(receiver_time, **kwargs) #extended_time, approx_t=True)
        #      for GW in GWs if GW.source_type in interp_type]

        if GWs[0].source_type in interp_type:
            hphc_call = 'interp_hphc'
            if ('precomputed' in kwargs and kwargs['precomputed']) or GWs[0].source_type == 'Numeric':
                pass
            else:
                jk = [GW.compute_hphc_td(receiver_time, **kwargs) for GW in GWs]
        else:
            hphc_call = 'compute_hphc_td'


        ### Compute GW effect on links
        self.t_min = t_min
        self.t_max = t_max
        self.dt = dt
        if 'tt' not in self.__dict__ or 'pos' not in self.__dict__:
            self.init_links(receiver_time, order=tt_order, **kwargs)

        nArms = self.orbits.number_of_arms if ilink is None else 1
        links = range(nArms) if ilink is None else [ilink]
        self.yArm = np.zeros([len(receiver_time), nArms])
        for link in links:
            for GW in GWs:
                self.yArm[:,link] += self._arm_response(receiver_time, link,
                                                        GW, hphc_call, **kwargs)
        return self.yArm


    def _arm_response(self, receiver_time, link, GW, call_hphc="compute_hphc", **kwargs):
        """Compute GW signal at receiver time for a single source and link.

        call_hphc can be switched to interpolation instead of actual
        hphc computation.

        """

        # Arm geometry
        receiver, emitter = self.orbits.get_pairs()[link]
        rem, rrec = self.pos[emitter-1,:,:], self.pos[receiver-1,:,:]
        # Choose between numpy or cupy
        xp = cp.get_array_module(rem) if cp is not None else np

        receiver_time = xp.array(receiver_time)
        tem = receiver_time - self.tt[:,link]
        r_ij = rrec - rem
        n = r_ij/xp.linalg.norm(r_ij, axis=0) # unit vector between receiver and emitter

        k,v,u = GW.basis
        k = xp.array(k)
        v = xp.array(v)
        u = xp.array(u)

        un, vn, kn = xp.dot(u, n), xp.dot(v, n), xp.dot(k, n)
        xi_p = 0.5*(un*un - vn*vn)
        xi_c = un*vn

        # Dot products of GW wave propagation vector and emitter, receiver postiional vectors
        kep = xp.dot(k, rem)
        krp = xp.dot(k, rrec)

        te = tem - kep
        tr = receiver_time - krp

        func = getattr(GW, call_hphc, None)
        if func is not None:
            # filter the kwargs belonging to func
            hpc_kwargs = {key: value for key, value in kwargs.items() if key in func.__code__.co_varnames}
            # Compute polarizations
            hpe, hce = func(te, **hpc_kwargs)
            hpr, hcr = func(tr, **hpc_kwargs)

        # Projected strain
        y = ((hpe - hpr)*xi_p + (hce - hcr)*xi_c )
        y/= (1.0-kn)

        # Move array from device(cupy) to host(numpy)
        if xp == cp:
            y = xp.asnumpy(y)

        return y


    def from_file(self, hdf5_filename):
        """ Load projected strain from file.
        """
        yArm, source_names, links, tmin, tmax, dt = from_file(hdf5_filename)
        self.yArm = yArm
        self.t_min = tmin
        self.t_max = tmax
        self.dt = dt
        self.source_names = source_names
        expected_links = ["%d-%d"%(a,b) for a,b in self.orbits.get_pairs()]
        assert links == expected_links

    def set_link(self):
        """ Define link names and indices.
        """
        self.links = ["%d-%d"%(a,b) for a,b in self.orbits.get_pairs()]
        self.dlink = dict(zip(self.links, range(len(self.links))))

    def to_file(self, hdf5_filename, fmt='ldc'):
        """Save projected strain on disk in hdf5 file.

         Args:
             fmt: string defining file format, either 'ldc' (default)
             or 'sim' to comply with LISA simulation tools.

        >>> proj = ProjectedStrain(orbits)
        >>> yArm = proj.arm_response(0, 1e4, 10, [GB])
        >>> proj.to_file("/tmp/test.hdf5")

        """
        if fmt == 'ldc':
            to_file(hdf5_filename, self.yArm, self.source_names, self.links,
                    self.t_min, self.t_max, self.dt)
        elif fmt == 'sim':
            with  h5py.File(hdf5_filename, 'w') as hdf5:
                t = np.arange(self.t_min, self.t_max, self.dt)
                hdf5.create_dataset('t', data=t)
                hdf5.attrs['gw_count'] = len(self.source_names)
                for j, link in enumerate(self.links):
                    dname = f'l_{link[0]}{link[-1]}'
                    hdf5.create_dataset(dname, data=self.yArm[:,j])
        else:
            print(f"Unkown format {fmt}")

    def slr2ilink(self, slr, perm=0):
        """Convert slr string (like 1-2) into link index between 0 and 6.

        >>> proj = ProjectedStrain(orbits)
        >>> proj.slr2ilink("132")
        2
        """
        fr = int(slr[0]) + perm
        to = int(slr[-1])+ perm
        fr = fr if fr<4 else fr%3
        to = to if to<4 else to%3
        key = "%d-%d"%(to, fr)
        return self.dlink[key]

    def interp(self, ilink):
        """ Interpolate projected strain.
        """
        return spline(np.arange(self.t_min, self.t_max, self.dt),
                      self.yArm[:,ilink])

    def y_slr(self, t, slr, d=[], perm=0, interp=None):
        """ Return TDI terms, following MLDC conventions (see MLDC user manual).

        y_{slr,d}(t) = y_slr (t - L_d) = y_slr (t - L_a - L_b ...)
        """
        # sum delay if any
        #delay = np.array([self.orbits.compute_travel_time(int(link[0]),
        #                                                  int(link[-1]), t)
        #                  for link in d]).sum(axis=0) if d else 0.
        delay = np.array([di for di in d]).sum(axis=0) if d else 0.
        t_delayed = t - delay

        # apply slr interpolator on delayed time range.
        if interp is None:
            ilink = self.slr2ilink(slr, perm=perm)
            interp = self.interp(ilink)
        res = interp(t_delayed)
        #res = self.sply[self.slr2ilink(slr)](t_delayed)
        return res

    def __permute(self, n):
        """ Return permuted indice n=0,1,2
        """

    def compute_tdi_x(self, t, tdi2=False, tt_order=0):
        """ Compute TDI X.
        """
        if tdi2:
            return self._compute_tdi_2(t, perm=0, tt_order=tt_order)
        return self._compute_tdi(t, perm=0, tt_order=tt_order)
    def compute_tdi_y(self, t, tdi2=False, tt_order=0):
        """ Compute TDI Y.
        """
        if tdi2:
            return self._compute_tdi_2(t, perm=1, tt_order=tt_order)
        return self._compute_tdi(t, perm=1, tt_order=tt_order)
    def compute_tdi_z(self, t, tdi2=False, tt_order=0):
        """ Compute TDI Z.
        """
        if tdi2:
            return self._compute_tdi_2(t, perm=2, tt_order=tt_order)
        return self._compute_tdi(t, perm=2, tt_order=tt_order)

    def _compute_tdi(self, t, perm=0, tt_order=0):
        """ Return X in time domain.

        We follow slr convention here.
        """
        delays = dict()
        for d in ["2-1", "1-3", "3-1", "1-2"]:
            fr = int(d[0]) + perm
            to = int(d[-1])+ perm
            fr = fr if fr<4 else fr%3
            to = to if to<4 else to%3
            delays[d] = self.orbits.compute_travel_time(fr, to, t, order=tt_order)

        spl = self.interp(self.slr2ilink("132", perm=perm))
        X = self.y_slr(t, "132", [delays["2-1"], delays["1-3"], delays["3-1"]], perm=perm, interp=spl)
        X-= self.y_slr(t, "132", [delays["2-1"]], perm=perm, interp=spl)
        spl = self.interp(self.slr2ilink("231", perm=perm))
        X+= self.y_slr(t, "231", [delays["1-3"], delays["3-1"]], perm=perm, interp=spl)
        X-= self.y_slr(t, "231", perm=perm, interp=spl)
        spl = self.interp(self.slr2ilink("123", perm=perm))
        X+= self.y_slr(t, "123", [delays["3-1"]], perm=perm, interp=spl)
        X-= self.y_slr(t, "123", [delays["3-1"], delays["1-2"], delays["2-1"]], perm=perm, interp=spl)
        spl = self.interp(self.slr2ilink("321", perm=perm))
        X+= self.y_slr(t, "321", perm=perm, interp=spl)
        X-= self.y_slr(t, "321", [delays["1-2"], delays["2-1"]], perm=perm, interp=spl)
        return X


    def _compute_tdi_2(self, t, perm=0, tt_order=0):
        delays = dict()
        for d in ["2-1", "1-3", "3-1", "1-2"]:
            fr = int(d[0]) + perm
            to = int(d[-1])+ perm
            fr = fr if fr<4 else fr%3
            to = to if to<4 else to%3
            delays[d] = self.orbits.compute_travel_time(fr, to, t, order=tt_order)

        spl = self.interp(self.slr2ilink("132", perm=perm))
        X = self.y_slr(t, "132", [delays["2-1"], delays["1-3"], delays["3-1"]],
                       perm=perm, interp=spl)
        X-= self.y_slr(t, "132", [delays["2-1"]], perm=perm, interp=spl)
        X+= self.y_slr(t, "132", [delays["3-1"], delays["1-3"], delays["2-1"],
                                  delays["1-2"], delays["2-1"]], perm=perm, interp=spl)
        X-= self.y_slr(t, "132", [delays["2-1"], delays["1-2"], delays["3-1"],
                                  delays["1-3"], delays["3-1"], delays["1-3"],
                                  delays["2-1"]], perm=perm, interp=spl)

        spl = self.interp(self.slr2ilink("231", perm=perm))
        X+= self.y_slr(t, "231", [delays["1-3"], delays["3-1"]], perm=perm, interp=spl)
        X-= self.y_slr(t, "231", perm=perm, interp=spl)
        X+= self.y_slr(t, "231", [delays["3-1"], delays["1-3"],
                                  delays["2-1"], delays["1-2"]], perm=perm, interp=spl)
        X-= self.y_slr(t, "231", [delays["2-1"], delays["1-2"], delays["3-1"],
                                  delays["1-3"], delays["3-1"], delays["1-3"]],
                       perm=perm, interp=spl)

        spl = self.interp(self.slr2ilink("123", perm=perm))
        X+= self.y_slr(t, "123", [delays["3-1"]], perm=perm, interp=spl)
        X-= self.y_slr(t, "123", [delays["3-1"], delays["1-2"], delays["2-1"]],
                       perm=perm, interp=spl)
        X+= self.y_slr(t, "123", [delays["3-1"], delays["1-3"], delays["2-1"],
                                  delays["1-2"], delays["2-1"], delays["1-2"],
                                  delays["3-1"]], perm=perm, interp=spl)
        X-= self.y_slr(t, "123", [delays["2-1"], delays["1-2"], delays["3-1"],
                                  delays["1-3"], delays["3-1"]], perm=perm, interp=spl)

        spl = self.interp(self.slr2ilink("321", perm=perm))
        X+= self.y_slr(t, "321", perm=perm, interp=spl)
        X-= self.y_slr(t, "321", [delays["1-2"], delays["2-1"]], perm=perm, interp=spl)
        X+= self.y_slr(t, "321", [delays["3-1"], delays["1-3"], delays["2-1"],
                                  delays["1-2"], delays["2-1"], delays["1-2"]],
                       perm=perm, interp=spl)
        X-= self.y_slr(t, "321", [delays["2-1"], delays["1-2"], delays["3-1"],
                                  delays["1-3"]], perm=perm, interp=spl)
        return X


if __name__ == "__main__":

    from astropy import units as un
    import doctest
    from ldc.waveform.waveform import HpHc
    from ldc.lisa.orbits import Orbits

    config = dict({"nominal_arm_length":2.5e9*un.m,
                   "initial_rotation":0*un.rad,
                   "initial_position":0*un.rad,
                   "orbit_type":"analytic"})

    pGB = dict({'Amplitude': 1.07345e-22,
                'EclipticLatitude': 0.312414*un.rad,
                'EclipticLongitude': -2.75291*un.rad,
                'Frequency': 0.00135962*un.Hz,
                'FrequencyDerivative': 8.94581279e-19*un.Unit('Hz2'),
                'Inclination': 0.523599*un.rad,
                'InitialPhase': 3.0581565*un.rad,
                'Polarization': 3.5621656*un.rad})

    GB = HpHc.type("my-galactic-binary", "GB", "TD_fdot")
    GB.set_param(pGB)
    orbits = Orbits.type(config)

    doctest.testmod()

#from lisanode import Graph
from lisanode.lisa.naming import indexed, fromto
from lisanode.lisa.naming import lisa_links, lisa_indices
from lisanode.lisa.config import ORBIT_TYPE, ORBIT_UPSAMPLING, ORBIT_SAMPLING_FREQ
from lisanode.compiler import Graph
from lisanode.lisa.config import LOW_SAMPLING_FREQ, HIGH_SAMPLING_FREQ, DOWNSAMPLING

#ORBIT_UPSAMPLING = 
ORBIT_SAMPLING_FREQ = 1/3600.#400.
ARM_LENGTH =  2.5E9
TT_ORDER = 2


class RefOrbits(Graph):

    def __init__(self):

        super().__init__("RefOrbits")

        self.add("KeplerianOrbits", name="orbits")

        # Set sampling parameters of orbits
        self.nodes["orbits"].fs = ORBIT_SAMPLING_FREQ
        self.nodes['orbits'].params['spacecraft_separation'] = ARM_LENGTH
        
        # Add and configure travel times
        self.add("TravelTimes", name="tt")
        if TT_ORDER==1:
            self.nodes["tt"].params = {
                'orderZero': True,
                'orderOneHalf': True,
                'orderOne': False,
            }
        elif TT_ORDER==2:
            self.nodes["tt"].params = {
                'orderZero': True,
                'orderOneHalf': True,
                'orderOne': True,
            }
        else:
            self.nodes["tt"].params = {
                'orderZero': True,
                'orderOneHalf': False,
                'orderOne': False,
            }


        
        # Connect orbits to travel times
        for i in range(3):
            self.connect("orbits.x[{}]".format(i), "tt.x[{}]".format(i))
            self.connect("orbits.y[{}]".format(i), "tt.y[{}]".format(i))
            self.connect("orbits.z[{}]".format(i), "tt.z[{}]".format(i))
            self.connect("orbits.vx[{}]".format(i), "tt.vx[{}]".format(i))
            self.connect("orbits.vy[{}]".format(i), "tt.vy[{}]".format(i))
            self.connect("orbits.vz[{}]".format(i), "tt.vz[{}]".format(i))
            # Publish orbits
            self.publish_output("orbits.x[{}]".format(i),"Orbitx[{}]".format(i))
            self.publish_output("orbits.y[{}]".format(i),"Orbity[{}]".format(i))
            self.publish_output("orbits.z[{}]".format(i),"Orbitz[{}]".format(i))
            self.publish_output("orbits.vx[{}]".format(i),"Orbitvx[{}]".format(i))
            self.publish_output("orbits.vy[{}]".format(i),"Orbitvy[{}]".format(i))
            self.publish_output("orbits.vz[{}]".format(i),"Orbitvz[{}]".format(i))
            
        #self.add("TravelTimes", name="travel_times")

        # Publish downsampled travel times
        for i, j in lisa_links():
            tt_name = fromto("tt", sending=i, receiving=j)
            travel_times = "tt.travel_times[{}][{}]".format(i-1, j-1)

            #self.add("Decimation<double>", "decimation_" + tt_name)
            #self.nodes["decimation_" + tt_name].downsampling = DOWNSAMPLING
            #self.connect(travel_times, "decimation_" + tt_name + ".input")
            #self.publish_output("decimation_" + tt_name + ".result", tt_name)
            self.publish_output(travel_times, tt_name)

class TestingLDCOrbits(Graph):

    def __init__(self):

        super().__init__("TestingLDCOrbits")


        #self.add("KeplerianOrbits", name="orbits")
        #self.nodes["orbits"].fs = ORBIT_SAMPLING_FREQ
        #self.nodes['orbits'].params['spacecraft_separation'] = ARM_LENGTH

        self.add("LDCOrbits", name="orbits")
        self.nodes["orbits"].fs = ORBIT_SAMPLING_FREQ
        self.publish_param("orbits" + ".arm_length", "arm_length")
        self.publish_param("orbits" + ".init_position", "init_position")
        self.publish_param("orbits" + ".init_rotation", "init_rotation")
        self.nodes['orbits'].params = {'init_position':0.0,
                                       'init_rotation':0.0, 'arm_length':ARM_LENGTH}
        
        
        # Add and configure travel times
        self.add("LDCTravelTimes", name="tt")
        #self.nodes["tt"].upsampling = ORBIT_UPSAMPLING
        self.publish_param("tt" + ".order", "order")
        self.nodes['tt'].params = {'order': TT_ORDER}
        
        # Connect orbits to travel times
        for i in range(3):
            self.connect("orbits.x[{}]".format(i), "tt.x[{}]".format(i))
            self.connect("orbits.y[{}]".format(i), "tt.y[{}]".format(i))
            self.connect("orbits.z[{}]".format(i), "tt.z[{}]".format(i))
            self.connect("orbits.vx[{}]".format(i), "tt.vx[{}]".format(i))
            self.connect("orbits.vy[{}]".format(i), "tt.vy[{}]".format(i))
            self.connect("orbits.vz[{}]".format(i), "tt.vz[{}]".format(i))
            # Publish orbits
            self.publish_output("orbits.x[{}]".format(i),"Orbitx[{}]".format(i))
            self.publish_output("orbits.y[{}]".format(i),"Orbity[{}]".format(i))
            self.publish_output("orbits.z[{}]".format(i),"Orbitz[{}]".format(i))
            self.publish_output("orbits.vx[{}]".format(i),"Orbitvx[{}]".format(i))
            self.publish_output("orbits.vy[{}]".format(i),"Orbitvy[{}]".format(i))
            self.publish_output("orbits.vz[{}]".format(i),"Orbitvz[{}]".format(i))
            
        # Publish downsampled travel times
        for i, j in lisa_links():
            tt_name = fromto("tt", sending=i, receiving=j)
            travel_times = "tt.travel_times[{}][{}]".format(i-1, j-1)
            self.publish_output(travel_times, tt_name)


class TestingTravelTimes_old(Graph):
    """ LISANode orbits and LDC travel times. 
    
    should give exact same results as RefOrbits. 
    """ 
    
    def __init__(self):

        super().__init__("TestingTravelTimes_old")


        self.add("KeplerianOrbits", name="orbits")
        self.nodes["orbits"].fs = ORBIT_SAMPLING_FREQ
        self.nodes['orbits'].params['spacecraft_separation'] = ARM_LENGTH

        # Add and configure travel times
        self.add("LDCTravelTimes", name="tt")
        #self.nodes["tt"].upsampling = ORBIT_UPSAMPLING
        self.publish_param("tt" + ".order", "order")
        self.nodes['tt'].params = {'order': TT_ORDER}
        
        # Connect orbits to travel times
        for i in range(3):
            self.connect("orbits.x[{}]".format(i), "tt.x[{}]".format(i))
            self.connect("orbits.y[{}]".format(i), "tt.y[{}]".format(i))
            self.connect("orbits.z[{}]".format(i), "tt.z[{}]".format(i))
            self.connect("orbits.vx[{}]".format(i), "tt.vx[{}]".format(i))
            self.connect("orbits.vy[{}]".format(i), "tt.vy[{}]".format(i))
            self.connect("orbits.vz[{}]".format(i), "tt.vz[{}]".format(i))
            # Publish orbits
            self.publish_output("orbits.x[{}]".format(i),"Orbitx[{}]".format(i))
            self.publish_output("orbits.y[{}]".format(i),"Orbity[{}]".format(i))
            self.publish_output("orbits.z[{}]".format(i),"Orbitz[{}]".format(i))
            self.publish_output("orbits.vx[{}]".format(i),"Orbitvx[{}]".format(i))
            self.publish_output("orbits.vy[{}]".format(i),"Orbitvy[{}]".format(i))
            self.publish_output("orbits.vz[{}]".format(i),"Orbitvz[{}]".format(i))
            
        # Publish downsampled travel times
        for i, j in lisa_links():
            tt_name = fromto("tt", sending=i, receiving=j)
            travel_times = "tt.travel_times[{}][{}]".format(i-1, j-1)
            self.publish_output(travel_times, tt_name)

class TestingTravelTimes(Graph):

    def __init__(self):

        super().__init__("TestingTravelTimes")


        #self.add("KeplerianOrbits", name="orbits")
        #self.nodes["orbits"].fs = ORBIT_SAMPLING_FREQ
        #self.nodes['orbits'].params['spacecraft_separation'] = ARM_LENGTH

        self.add("LDCOrbits", name="orbits")
        self.nodes["orbits"].fs = ORBIT_SAMPLING_FREQ
        self.publish_param("orbits" + ".arm_length", "arm_length")
        self.publish_param("orbits" + ".init_position", "init_position")
        self.publish_param("orbits" + ".init_rotation", "init_rotation")
        self.nodes['orbits'].params = {'init_position':0.0,
                                       'init_rotation':0.0, 'arm_length':ARM_LENGTH}

        # Add and configure travel times
        self.add("TravelTimes", name="tt")
        if TT_ORDER==1:
            self.nodes["tt"].params = {
                'orderZero': True,
                'orderOneHalf': True,
                'orderOne': False,
            }
        elif TT_ORDER==2:
            self.nodes["tt"].params = {
                'orderZero': True,
                'orderOneHalf': True,
                'orderOne': True,
            }
        else:
            self.nodes["tt"].params = {
                'orderZero': True,
                'orderOneHalf': False,
                'orderOne': False,
            }

        #self.nodes["tt"].upsampling = ORBIT_UPSAMPLING
        #self.publish_param("tt" + ".order", "order")
        #self.nodes['tt'].params = {'order': TT_ORDER}
        
        # Connect orbits to travel times
        for i in range(3):
            self.connect("orbits.x[{}]".format(i), "tt.x[{}]".format(i))
            self.connect("orbits.y[{}]".format(i), "tt.y[{}]".format(i))
            self.connect("orbits.z[{}]".format(i), "tt.z[{}]".format(i))
            self.connect("orbits.vx[{}]".format(i), "tt.vx[{}]".format(i))
            self.connect("orbits.vy[{}]".format(i), "tt.vy[{}]".format(i))
            self.connect("orbits.vz[{}]".format(i), "tt.vz[{}]".format(i))
            # Publish orbits
            self.publish_output("orbits.x[{}]".format(i),"Orbitx[{}]".format(i))
            self.publish_output("orbits.y[{}]".format(i),"Orbity[{}]".format(i))
            self.publish_output("orbits.z[{}]".format(i),"Orbitz[{}]".format(i))
            self.publish_output("orbits.vx[{}]".format(i),"Orbitvx[{}]".format(i))
            self.publish_output("orbits.vy[{}]".format(i),"Orbitvy[{}]".format(i))
            self.publish_output("orbits.vz[{}]".format(i),"Orbitvz[{}]".format(i))
            
        # Publish downsampled travel times
        for i, j in lisa_links():
            tt_name = fromto("tt", sending=i, receiving=j)
            travel_times = "tt.travel_times[{}][{}]".format(i-1, j-1)
            self.publish_output(travel_times, tt_name)

from lisanode.compiler import Graph

class Response(Graph):
    def __init__(self):
        
        super().__init__("Response")
        LISA_PHYSICS_FS = 3. # Hz
        
        self.add("ReadHDF5<double>", name="h5")
        self.nodes["h5"].params = {'path': "input.h5",
                                   'dataset_path': "strain"}
        self.nodes["h5"].fs = LISA_PHYSICS_FS
        
        avg =  "avg"
        self.add("EllipticFilter<double>", avg)
        
        # target_dt = 5 
        # self.nodes[avg].params = {
        #     'passband_freq': 0.4 * (1/target_dt), 
        #     'stopband_freq': 0.4985 * (1/target_dt), 
        #     'minimum_passband_gain': 0.1,
        #     'minimum_attenuation': 100}
        target_dt = 5 
        self.nodes[avg].params = {
            'passband_freq': 0.40 * (1/3.45), #3.65
            'stopband_freq': 0.48 * (1/3.45), #3.65
            'minimum_passband_gain': 0.1,
            'minimum_attenuation': 100}

        
        decim = "decimation" 
        self.add("Decimation<double>", decim)
        self.nodes[decim].downsampling =int(LISA_PHYSICS_FS/(1./target_dt))
        self.nodes[decim].lag = 10
        self.connect("h5.data[0]", avg + ".input")
        self.connect(avg+".result", decim+".input")
        self.publish_output(decim + ".result", "out")
        self.publish_output(avg + ".result", "raw_out")

# nfc2
new repository for nfc

This is a public repository for Numeric Fluid Channel - a distributed parallel LBM-based code 
for CFD calculations on a "channel" geometry: parallelepiped fluid channel with constant velocity
inlet boundary condition and constant pressure outlet boundary condition.  

Acceleration is provided via GPUs with OpenACC directives.  Non-GPU equipped nodes benefit from
OpenMP with some modification (or just use the OpenMP branch.)

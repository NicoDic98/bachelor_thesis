import h5py

f = h5py.File('out.h5', 'r')

print(list(f.keys())[0])
level0 = f.get(list(f.keys())[0])
print(level0.keys())
meas = level0.get(list(level0.keys())[1])
print(meas.keys())
magn = meas.get(list(meas.keys())[0])
print(magn.attrs.keys())
print(dict(magn.attrs.items()))

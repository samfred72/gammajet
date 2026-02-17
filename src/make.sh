rm drawer.o ana.o histmaker.o object.o pho_object.o jet_object.o libgammajet.so
rm /home/samson72/root/lib/libgammajet.so
$(root-config --cxx) -c -fPIC drawer.cc ana.cc histmaker.cc object.cc pho_object.cc jet_object.cc `root-config --cflags`
echo ".o files made"
$(root-config --cxx) -shared -o libgammajet.so drawer.o ana.o histmaker.o object.o pho_object.o jet_object.o `root-config --libs`
echo "library made"
cp libgammajet.so /home/samson72/root/lib/

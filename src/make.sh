rm unfolder.o drawer.o ana.o histmaker.o object.o pho_object.o jet_object.o treeuser.o libgammajet.so
rm /home/samson72/root/lib/libgammajet.so
$(root-config --cxx) -c -fPIC -Wno-deprecated-declarations\
  drawer.cc \
  unfolder.cc \
  ana.cc \
  histmaker.cc \
  object.cc \
  pho_object.cc \
  jet_object.cc \
  treeuser.cc \
  `root-config --cflags`
echo ".o files made"
$(root-config --cxx) -shared -Wno-deprecated-declarations -o \
  libgammajet.so \
  drawer.o \
  unfolder.o \
  ana.o \
  histmaker.o \
  object.o \
  pho_object.o \
  jet_object.o \
  treeuser.o \
  `root-config --libs`
echo "library made"
cp libgammajet.so /home/samson72/root/lib/

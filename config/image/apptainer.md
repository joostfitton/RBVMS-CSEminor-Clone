# Apptainer

# Run an image
apptainer shell rbvms.sif

# Use an overlay
apptainer overlay create --fakeroot --size 1024 overlay.img
apptainer shell --fakeroot --overlay overlay.img rbvm.sif

# Build an image
apptainer build rbvms.sif rbvms.def




# Run a sanbox image with persistent changes
apptainer shell --writable rbvms/


# Build a sandbox image
apptainer build --sandbox rbvms rbvms.def




# Install apptainer
https://apptainer.org/docs/admin/main/installation.html#install-ubuntu-packages

sudo apt update
sudo apt install -y software-properties-common
sudo add-apt-repository -y ppa:apptainer/ppa
sudo apt update
sudo apt install -y apptainer
sudo apt install -y apptainer-suid



## Slambook2

Install dependences

	sudo apt-get install -y libfmt-dev

Udpade submodules (3rd party)

	git submodule update --init --recursive
	
Build submodules, examples Sophus

	cd Sophus
	git checkout 785fef35b7d9e0fc67b4964a69124277b7434a44   (master ok ,    13fb3288311485dc94e3226b69c9b59cd06ff94e is old ko   )
	
	NO NEEDS TO COMPILE
	mkdir build && cd build && cmake .. -j$(nproc) && make    (OPTIONAL)

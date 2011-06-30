all:
	make -C ../doc html
	scp -r ../doc/_build/html/* spin:Sites/rmt

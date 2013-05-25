PERF_MAX_LOG_N=8

OPTIONS=-Iinclude

%.o : src/%.c
	gcc $(OPTIONS) -o $@ -c $<

bin/quantify : src/quantify.c geneindex.o
	gcc $(OPTIONS) -o $@ $^

bin/testgeneindex : src/testgeneindex.c geneindex.o
	gcc $(OPTIONS) -o $@ $^

cover : 
	nosetests --with-coverage --cover-html --cover-package pade

test : bin/testgeneindex
	bin/testgeneindex
#	nosetests --with-doctest

clean :
	rm -f *.log tests/*~ tests/*.pyc site.tar
	cd doc; make clean
	rm -rf doc/html/generated
	rm -rf cover
	rm -f `find . -name \*~`
	rm -f `find . -name \*.pyc`

site :
	rm -rf doc/generated
	sphinx-apidoc pade -o doc/generated
	cd doc; make clean html

deploy_site:

	cd doc; make html
	cd doc/_build/html; tar cf ../../../site.tar *
	git checkout gh-pages
	tar xf site.tar
	git add `tar tf site.tar`
	git commit -m 'Updating site'
	git push origin gh-pages
	git checkout master


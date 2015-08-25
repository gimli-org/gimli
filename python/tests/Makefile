all: test

tests: test

clean:
	rm -rf *.vector *.bmat *.bms *.matrix __pycache__

test:
	python -c "import pygimli; pygimli.test(show=False)"
	py.test

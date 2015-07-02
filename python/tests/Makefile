all: test

clean:
	rm -rf *.vector *.bmat *.bms *.matrix __pycache__

test:
	echo "####################"
	echo "Running all tests..."
	echo "####################"
	sleep 1
	python -c "import pygimli; pygimli.test(show=True)"
	py.test .

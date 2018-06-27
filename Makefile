coverage = coverage run $(omit) --source gffparser -m pytest

unit:
	pytest --durations=3 gffparser 

clean:
	find . -name '*.pyc' | xargs rm -f
	find . -name '__pycache__' | xargs rm -rf
	rm -rf cover .coverage
	rm -rf gffparser.egg-info

install:
	pip install -e .[testing]

coverage:
	rm -rf cover .coverage
	coverage run --source gffparser -m pytest gffparser
	coverage html -d cover
	coverage report 


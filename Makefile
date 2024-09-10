set_env:
	pyenv virtualenv 3.10.6 ddapp
	pyenv local ddapp

reinstall_package:
	@pip uninstall -y ddapp || :
	@pip install -e .

reinstall_packagex:
	@brew install pipx
	@pipx uninstall -y ddapp || :
	@pipx install -e .

streamlit:
	-@streamlit run app.py

clean:
	@rm -fr */__pycache__
	@rm -fr __init__.py
	@rm -fr build
	@rm -fr dist
	@rm -fr *.dist-info
	@rm -fr *.egg-info
	-@rm model.joblib

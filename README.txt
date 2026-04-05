To set up environment:

While you're in the folder with your project run (in this exact order):

1. Set-ExecutionPolicy -ExecutionPolicy RemoteSigned -Scope CurrentUser -- tells windows to mind its own business
2. python -m venv .venv -- creates the virtual environment where all of your packages will be installed into
3. .venv\Scripts\Activate.ps1 -- activates the virtual environment, connects ur python and terminal basically
4. pip install -r requirements.txt -- installs all of the packages listed in the txt file to your virtual environment so python understands your code
 
Happy Coding!
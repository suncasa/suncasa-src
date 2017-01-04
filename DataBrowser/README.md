# DataBrowser
A interactive tool for solar radio imaging spectroscopy data analysis

#Build Requirements
```
bokeh>=0.12.0
sunpy>=0.7.4
pandas>=0.19.2
scikit-image>=0.12.3
```

# Usage
Install [dependencies](requirements.txt) using:

```
$ pip install -r requirements.txt
```

set the PATH environment variable.
If in C-Shell, set your ~/.cshrc by adding:

```swift
setenv SUNCASA "<your suncasa path>"
setenv SUNCASADB "<your database path>"
setenv SUNCASAPY46 "<your CASA 4.6 path>"
setenv SUNCASAPY47 "<your CASA 4.7 path>"
alias browser csh $SUNCASA/DataBrowser/browser_csh.sh
```

If in B-Shell, set your ~/.bashrc by adding:
```swift
export SUNCASA = "<your suncasa path>"
export SUNCASADB = "<your database path>"
export SUNCASAPY46 = "<your CASA 4.6 path>"
export SUNCASAPY47 = "<your CASA 4.7 path>"
alias browser='bash $SUNCASA/DataBrowser/browser_bsh.sh'
```

To launch the Databrowser, type in browser in your terminal.

# Authors
### Original Author and Development Lead
- Sijie Yu ([@sjyu1988](https://github.com/sjyu1988))

### Co-Author
- Bin Chen([@binchensun](https://github.com/binchensun))


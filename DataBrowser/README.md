# DataBrowser
A interactive tool for solar radio imaging spectroscopy data analysis

# Usage
set the PATH environment variable.
If in C-Shell, set your ~/.cshrc by adding:

```swift
setenv SUNCASA "<your suncasa path>"
setenv SUNCASADB "<your database path>"
setenv SUNCASAPY46 "<your CASA 4.6 path>"
setenv SUNCASAPY47 "<your CASA 4.7 path>"
alias browser bokeh serve --show $SUNCASA/DataBrowser/EvtBrowser --port `(netstat  -atn | awk '{printf "%s\n%s\n", $4, $4}' | grep -oE '[0-9]*$'; seq 32768 61000) | sort -n | uniq -u | head -n 1`
```

If in B-Shell, set your ~/.bashrc by adding:
```swift
export SUNCASA = "<your suncasa path>"
export SUNCASADB = "<your database path>"
export SUNCASAPY46 = "<your CASA 4.6 path>"
export SUNCASAPY47 = "<your CASA 4.7 path>"
alias browser='bokeh serve --show $SUNCASA/DataBrowser/EvtBrowser --port `port=32768; while netstat -atn | grep -q $port; do port=$(expr $port + 1); done; echo $port`'
```

To launch the Databrowser, type in browser in your terminal.

# Authors
### Original Author and Development Lead
- Sijie Yu ([@sjyu1988](https://github.com/sjyu1988))

### Co-Author
- Bin Chen([@binchensun](https://github.com/binchensun))

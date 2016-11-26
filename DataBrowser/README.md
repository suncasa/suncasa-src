# DataBrowser
A interactive tool for radio imaging spectroscopy data analysis

# Usage
set the PATH environment variable.
If in C-Shell, set your ~/.cshrc by adding:

```swift
setenv SUNCASA "<your suncasa path>"
setenv SUNCASADB "<your database path>"
setenv SUNCASAPY46 "<your CASA 4.6 path>"
setenv SUNCASAPY47 "<your CASA 4.7 path>"
alias browser bokeh serve --show $SUNCASA/DataBrowser/EvtBrowser --port `netstat -atn | awk ' /tcp/ {printf("%s\\\\\\\\\\\\\\\\n",substr($4,index($4,":")+1,length($4) )) }' | sed -e "s/://g" | sort -rnu | awk '{array [$1] = $1} END {i=32768; again=1; while (again == 1) {if (array[i] == i) {i=i+1} else {print i; again=0}}}'`
```

If in B-Shell, set your ~/.bashrc by adding:
```swift
export SUNCASA = "<your suncasa path>"
export SUNCASADB = "<your database path>"
export SUNCASAPY46 = "<your CASA 4.6 path>"
export SUNCASAPY47 = "<your CASA 4.7 path>"
alias browser='bokeh serve --show $SUNCASA/DataBrowser/EvtBrowser'
```

To launch the Databrowser, type in browser in your terminal.

# Authors
### Original Author and Development Lead
- Sijie Yu ([@sjyu1988](https://github.com/sjyu1988))

### Co-Author
- Bin Chen([@binchensun](https://github.com/binchensun))

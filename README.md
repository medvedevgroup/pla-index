# pla-index with strobealign
We replace `randstrobe_start_indices` with `repeat-stretched pla-index`.
To run strobealign with pla-index, CMakeLists.txt needs to be updated with include and library folder of SDSL.\
To update CMakeLists.txt run the following from the current directory:

```
./update_config.sh SDSL_INCLUDE_PATH SDSL_LIB_PATH

Parameter description:
SDSL_INCLUDE_PATH: path of the SDSL include folder
SDSL_LIB_PATH: path of the SDSL library folder
```

Afterwards, strobealign can be run normally with `--eps [INT]` parameter to specify the epsilon value to be used in the pla-index.
Details are specified in the README file inside the strobealign folder.
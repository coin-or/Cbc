# script to format source code and include (if not included yet),
# vim formatting settings (which are automatically loaded by vim)
# Haroldo - 2019

for file in *.[ch]pp *.h;
do
    sourceName=`basename $file`
    echo formatting "$sourceName"
    clang-format -i -style=file $file

    # adding vim modeline if not included yet
    if ! grep -q "/* vi: softtabstop=" $file; then
      echo '' >> $file
      echo  '/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2' >> $file
      echo  '*/' >> $file
    fi
done

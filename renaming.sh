newName=$1
echo "renaming to "$newName
git mv shablona $newName
git mv $newName/shablona.py $newName/$newName.py
git mv $newName/tests/test_shablona.py $newName/tests/test_$newName.py

git commit -a -m "Moved names from 'shablona' to "$newName

git rm $newName/data/*
git commit -a -m "Removed example 'shablona' data"

git grep -l 'shablona' | xargs sed -i 's/shablona/'$newName'/g'

git rm doc/conf.py

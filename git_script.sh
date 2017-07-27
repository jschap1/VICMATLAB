#! /bin/bash

# This is a script to push changes to the online git repository.

git init
git add .
git commit -m "Update"
git push -u origin master
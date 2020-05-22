#! /bin/bash

# configure your name and email if you have not done so
git config --global user.email "benkeser@emory.edu"
git config --global user.name "David Benkeser"

# clone the repository to the book-output directory
echo $TRAVIS_REPO_SLUG
git clone -b gh-pages \
  https://${GITHUB_PAT}@github.com/${TRAVIS_REPO_SLUG}.git \
  slapnap
cd slapnap
ls 
git rm -rf *
cp -r ../docs/_book/* ./
git add --all *
git commit -m "Update the book"
git push -q origin gh-pages
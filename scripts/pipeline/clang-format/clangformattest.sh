#!/bin/bash

## Javier Galan got this script from
## https://gitlab.cern.ch/VecGeom/VecGeom/commit/8c723ded243dbfaac59d372b524366db4577a728

#
# This is running a clang-format test
# by doing a filtering step and then analysing
# the result of applying ./scripts/clang-format-and-fix-macros.sh
#

# check that we are in a clean state in order to prevent accidential
# changes
cleanstate=`git status | grep "modified"`
if ! [[ -z $cleanstate ]]; then
	echo "Script must be applied on a clean git state"
	return 0
fi


cd ../../../
find . \( -name '*.h' -or -name '*.cxx' -or -name '*.cc' -or -name '*.C' \) -print0 | xargs -0 clang-format -i

# check if something was modified
notcorrectlist=`git status | grep "modified"`
if [[ -z $notcorrectlist ]]; then
	echo "Excellent. **VERY GOOD FORMATTING!** :thumbsup:"
	exit 0;
else
	echo "The following files have clang-format problems (showing patches)";
	for f in $notcorrectlist; do
		echo $f
		#git diff $f
	done

	#We now repair commiting
	git config --global user.email "runner@lfna.unizar.es"
	git config --global user.name "runner"
	git config --global push.default simple
	git remote set-url --push origin https://runner:uwh7Ui*087@lfna.unizar.es/iaxo/RestAxionLib.git
	git add -u
	git commit -m "Pipeline clang-format automatic execution"
	git status
	git config --global http.sslverify "false"
    echo "git push origin HEAD:$CI_COMMIT_REF_NAME"
	git push origin HEAD:$CI_COMMIT_REF_NAME
	echo "Clang-format should have generated a commit to fix code formatting"
	echo "This validation should be solved in the next run."
	exit 0;
fi

# cleanup changes in git
git reset HEAD --hard

exit 1


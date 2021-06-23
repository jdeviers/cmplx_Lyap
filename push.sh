echo "call like: ./push.sh \"COMMIT MESSAGE\", WITH DOUBLE QUOTES."
[[ "$1" == "" ]] && echo -e "\n*** ERROR ***\nEmpty message: quitting." && exit

git config --global user.email "jd756@exeter.ac.uk"
git config --global user.name "jdeviers_SSD"

git add .
git commit -m "${message}"

git push origin master

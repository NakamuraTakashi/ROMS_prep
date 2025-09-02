#!/bin/sh
#
# 1. Please visit CDS web page (https://cds.climate.copernicus.eu/)
# 2. Login and visit:
#     https://cds.climate.copernicus.eu/how-to-api
# 3. See "Linux users" instruction, and check your "key: ", 
#    and replace "********-****-****-****-************" 
#    with your key.
# 4. Rename this "setup.sh" to "my_setup.sh".
#     *Should not upload the "my_setup.sh" file to any open repositories, such as GitHub,
#      for keeping the secret of your key! 
# 5. Change the permission of this file to -rwx------ (700).
#    > chmod 700 my_setup.sh
# 6. Run the "my_setup.sh" as follows:
#    > ./my_setup.sh
#
echo "url: https://cds.climate.copernicus.eu/api" > $HOME/.cdsapirc
echo "key: ********-****-****-****-************" >> $HOME/.cdsapirc
#
#pip3 install cdsapi
pip install "cdsapi>=0.7.2" --break-system-packages

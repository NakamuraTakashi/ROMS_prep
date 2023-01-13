#!/bin/sh
#
# 1. Please visit CDS web page (https://cds.climate.copernicus.eu/cdsapp#!/home)
# 2. Login and visit:
#     https://cds.climate.copernicus.eu/api-how-to#use-the-cds-api-client-for-data-access
# 3. Check your "key: ", and replace "******:************************************" 
#    with your key.
# 4. Rename this "setup.sh" to "my_setup.sh".
#     *Should not upload the "my_setup.sh" file to any open repositories, such as GitHub,
#      for keeping the secret of your key! 
# 5. Run the "my_setup.sh" as follows:
#    > ./my_setup.sh
#
echo "url: https://cds.climate.copernicus.eu/api/v2" > $HOME/.cdsapirc
echo "key: ******:************************************" >> $HOME/.cdsapirc
#
pip3 install cdsapi

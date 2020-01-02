# Matching genetic sampling locations to RAM Legacy fishery stocks
version 4/8/2016

We are trying to add a “stockid” to the genetic sampling locations in the microsatellite and mtDNA database that we have assembled. The stockid comes from the RAM Legacy database (http://www.ramlegacy.org), which records that status of many fish populations around the world (but not all). 

We have three files:
1. msat_to_match: the microsatellite data with an extra stockid column for us to fill
2. mtdna_to_match: the mtDNA data with an extra stockid column for us to fill
3. stocks_to_match: the table of RAM Legacy stocks

## Instructions
1. For each Site line in the msat_to_match file, try to find a matching stocklong in the stocks_to_match file for the same species (use the fbsci column to match species).

   Example: 
Brevoortia patronus is the first species (fbsci column) in the msat_to_match file, and the first location is “Galveston, TX” in “USA” (Site and County columns). We turn to the stocks_to_match file and find the lines for Brevoortia patronus. One of those lines (actually, the only line) says “Gulf menhaden Gulf of Mexico” in the stocklong column. Galveston, TX is in the Gulf of Mexico, so this is a match. 

2. If a match is found, copy the stockid from stocks_to_match to the appropriate line of the msat_to_match file.

Example continued:
Copy “MENATGM” from the stockid column of stocks_to_match into the first line of the stockid column of msat_to_match (the line for Brevoortia patronus, Galveston, TX).

3. If no match is found, type NA in the stockid column of the msat_to_match file.

4. Repeat for the mtdna_to_match file.

## Understanding stocklong names and matching locations
The stocklong column in stocks_to_match contains the common name for the species and a text description of the location of that population. Some locations are easy to understand, like Gulf of Mexico. Others are harder, like “NAFO 4X” or “ICES 22-24-IIIa.” The latter refer specifically to NAFO or ICES zones, and you can refer to maps of these zones to determine whether a particular location is inside those zones. The latitude and longitude from the msat_to_match and the mtdna_to_match files may also be helpful in understanding where a given Site is found. If any locations or stocklong’s are hard to understand, please email me.

## Logbook
Keeping a data entry logbook is very important. Each day of work, please make a note of the date, the hours worked, and the progress made. Most importantly, please write down any locations, sites, or matches that were confusing or that you are uncertain about. This allows us to go back and double check them together. Please also note any relevant thoughts, concerns, or particular methods that you use.

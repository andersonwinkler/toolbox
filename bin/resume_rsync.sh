#!/bin/sh
# Reliable file transfer
# Original script by Kleines Tourenbuch
#
# Added trap, parsable arguments and better feedback by Anderson Winkler
# Obs: to allow passwordless SSH login, use public keys (RSA)
# 09.Feb.2010

FROM=$1
DEST=$2

# Define a function for Ctrl+C
trap bashtrap INT
bashtrap()
{
  break
  exit 1
}

# Try rsync for $MAX_RESTARTS times 
I=0
MAX_RESTARTS=120000
LAST_EXIT_CODE=1
while [ ${I} -le ${MAX_RESTARTS} ]
do
  I=$(( ${I} + 1 ))
  echo "===[ Start of rsync: #${I} ]==="
  rsync -av --partial --progress --copy-links --timeout=120 -e "ssh" ${FROM} ${DEST}
  LAST_EXIT_CODE=$?
  if [ ${LAST_EXIT_CODE} -eq 0 ]; then
    break
  fi
  sleep 5
done

# Check if successful
if [ ${LAST_EXIT_CODE} -ne 0 ]; then
  echo "rsync failed for ${I} times. Giving up."
else
  echo "rsync successful after ${I} times."
fi



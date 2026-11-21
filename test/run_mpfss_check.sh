#!/bin/sh
# usage: run_mpfss_check.sh <binary> <port> <mode> <field> [<expected party 1|2> <expected message>]
# honest: both parties exit 0 and print "round passed".
# otherwise: the expected party exits non-zero and prints the message, and the
# other party does not complete the round as the receiver (BOB never passes).
B=$1; PORT=$2; MODE=$3; FIELD=$4; WHO=$5; MSG=$6
D=$(mktemp -d)
"$B" 1 "$PORT" "$MODE" "$FIELD" > "$D/1.log" 2>&1 & P1=$!
"$B" 2 "$PORT" "$MODE" "$FIELD" > "$D/2.log" 2>&1; R2=$?
wait $P1; R1=$?
echo "--- ALICE (exit $R1)"; cat "$D/1.log"; echo "--- BOB (exit $R2)"; cat "$D/2.log"
OK=1
if [ "$MODE" = honest ]; then
  [ $R1 -eq 0 ] && [ $R2 -eq 0 ] && grep -q "round passed" "$D/1.log" && grep -q "round passed" "$D/2.log" || OK=0
else
  eval R=\$R$WHO
  [ $R -ne 0 ] && grep -q "$MSG" "$D/$WHO.log" || OK=0
  grep -q "round passed" "$D/2.log" && OK=0   # the receiver must never accept a cheating round
fi
rm -rf "$D"
[ $OK -eq 1 ] && echo "RESULT: as expected" && exit 0
echo "RESULT: UNEXPECTED"; exit 1

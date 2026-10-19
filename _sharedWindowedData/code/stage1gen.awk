# Generalised stage-1 filter. Lazy parsing: never reference $1..$N so gawk does not
# split the multi-kilobyte blob fields on every row.
#   -v nskip=1  trained : Experiment,Trial,Cell,[Trace],TraceSpikes,...
#   -v nskip=2  naive   : Experiment,Trial,Cell,[layer,Trace],TraceSpikes,...
# Output: Experiment,Trial,Cell,Behav_Cond,Visual_Stim,TraceSpikes(hex)
BEGIN{
  if (nskip=="") nskip=1;
  while ((getline line < qcf) > 0) qc[line]=1;
  close(qcf);
  while ((getline line < ttf) > 0) { split(line, a, "\t"); tt[a[1]] = a[2] "," a[3]; }
  close(ttf);
  kept=0; n=0;
}
{
  n++;
  if (match($0, /^[^,]*,[^,]*,[^,]*,/) == 0) next;
  head = substr($0, 1, RSTART+RLENGTH-2);          # Experiment,Trial,Cell
  split(head, h, ",");
  if (!((h[1] "|" h[3]) in qc)) next;
  k2 = h[1] "|" h[2];
  if (!(k2 in tt)) next;
  rest = substr($0, RSTART+RLENGTH);               # fields after Cell
  for (i = 0; i < nskip; i++) {                    # drop layer / Trace
    p = index(rest, ",");
    if (p == 0) { rest = ""; break; }
    rest = substr(rest, p+1);
  }
  q = index(rest, ",");
  if (q > 0) ts = substr(rest, 1, q-1); else ts = rest;
  if (ts == "" || substr(ts,1,2) != "0x") next;    # skip NULL / malformed blobs
  kept++;
  print h[1] "," h[2] "," h[3] "," tt[k2] "," ts;
}
END { printf "[stage1gen] read %d rows, kept %d\n", n, kept > "/dev/stderr" }

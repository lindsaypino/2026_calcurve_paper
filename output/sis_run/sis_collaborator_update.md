**Subject:** 240610 SIS response curves — LOD/LOQ came out clean, but we never found the ceiling

Hi [Name],

Thanks again for sending over the 240610 SIS peptide-response data. I finally got it
through our LOD/LOQ pipeline this week, and I wanted to share what we found while it's
fresh — plus one ask at the end.

Short version: the curves are lovely and the detection and quantification limits came
out clean. But the dataset doesn't have an upper limit of quantitation hiding in it,
because none of the peptides ever actually saturate. To see the ceiling we'd need to
spike higher.

**The good news.** Once I pointed the algorithm at the heavy channel (light was all
#N/A, which is what you'd expect for the spiked standards) and handed it your
per-peptide concentration multipliers, it fit 19 of the 26 peptides with a finite LOD
and 18 with an LOQ. One thing that made me trust the output: it independently
reproduced your own QC notes. Every peptide you flagged in group_A as "no signal" came
back with an infinite LOD, and in group_B the only one it could fit was EASGLSADSLAR —
exactly the one you'd marked "ok". Two completely different routes to the same call,
which is about the most reassuring thing you can ask for. Your annotations saved me a
lot of squinting, so thank you for those.

(One inside-baseball note in case it's useful: I had to run with our noise-point
minimum set to zero. The high-multiplier groups are so far above the noise floor that
they never give the fit a noise plateau to anchor on — at our default setting only 6 of
26 peptides got an LOD. Not a problem with the data, just a knob.)

**The ULOQ problem, which is the ask.** Our ULOQ — upper limit of quantitation, i.e.
the concentration where the detector stops responding proportionally and the curve
flattens — needs to actually *see* that flat part, specifically at least two curve
points sitting up on the plateau. Right now zero of 26 peptides have that.

Here's the tell. With a 3-fold dilution, each step up the curve should multiply signal
by about 3x. Two points that have both hit the ceiling would come in around 1x. The
measured top steps are 2.2–2.8x. So there *is* a little compression at the top, which
is real and interesting in its own right, but the curves are still climbing hard at the
highest level you spiked. It's a bit like trying to find where a speaker starts to
distort when the dial only goes to 5 — nothing wrong with the measurement, we just
never got loud enough.

The one exception is LGQHLATEPLGTNSWER, which does visibly bend at the top (1.36x).
I've attached the figure — all 26 curves, log-log, one column per spike group, dashed
grey is ideal proportional response. That one peptide is the fourth panel down in the
group E column, and I think you'll spot it immediately.

So the ask: is it feasible to extend the top of the curve? My guess is two or three
more 3-fold steps above the current top would be enough to bracket the ceiling, and I'd
start with the group_E (x10,000) peptides since they're already closest. Tighter
spacing up top instead of wider range would also work, and might even be better — but
that's a real instrument-time question and you'd know better than I would whether it's
worth it.

Either way this wasn't wasted. "The Ultra's linear range extends past the top of a
5-point SIS curve" is a genuinely useful result and may well end up in the paper as a
supplement. I just don't want to build the ULOQ argument on a plateau we can't see.

Happy to jump on a call and walk through the figure if that's easier than email.

Thanks again,
Lindsay

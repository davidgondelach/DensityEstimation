% Test J2000 to TEME transformation

rr_J2000_test = [3961.7442603 6010.2156109 4619.3625758]';
vv_J2000_test = [-5.314643386 3.964357585 1.752939153]';

% June 28, 2000, 15: 8:51.655 000 UTC
year = 2000;
mon = 6;
day = 28;
hr = 15;
min = 8;
sec = 51.655;
timezone = 0;
[jdate, jdatefrac] = jday(2000, 6, 28, 15, 8, 51.655);

global EOPMat;
[ ~, ~, dut1, ~, ddpsi, ddeps, dat ] = computeEOP_Celestrak( EOPMat, jdate+jdatefrac );

[ut1, tut1, jdut1, jdut1frac, utc, tai, tt, ttt, jdtt,jdttfrac, tdb, ttdb, jdtdb,jdtdbfrac, tcg, jdtcg,jdtcgfrac, tcb, jdtcb,jdtcbfrac ] ...
    = convtime ( year, mon, day, hr, min, sec, timezone, dut1, dat );

[rteme1, vteme1, ateme1] = eci2teme(rr_J2000_test, vv_J2000_test, [0 0 0]', ttt, ddpsi, ddeps);

[rj2000, vj2000] = convertTEME2ECI(rteme1, vteme1, jdate+jdatefrac );


[rj2000v2, vj2000v2] = convertTEMEtoJ2000(rteme1, vteme1, jdate+jdatefrac);
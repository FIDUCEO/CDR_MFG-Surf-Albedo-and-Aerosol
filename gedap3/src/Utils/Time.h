/**********************************************************************
* The copyrights for the GEDAP algorithm and computer codes, remain with 
* Rayference SPRL as an Intellectual Property Right 
*
* Any use, in source and binary forms, distribution, modification, for 
* commercial, business, research and any other purposes is
*
*                           *** STRICTLY FORBIDDEN ***
*
* GEDAP may not be distributed or sold to any other commercial, 
* business, research or other partners under any circumstances.
* In any case this comment is part of the code as Software Legal Information,
* it cannot be modified and must be kept as header for all source codes in 
* which it appears.
*
* All documents and software are protected under the international Copyright 
* Laws, and Rayference SPRL reserves all rights.
*
*
* Author       :     Rayference Copyright (c) 
*
*********************************************************************/
#ifndef TILE_MAKER_MODULE_UTILS_TIME_H_
#define TILE_MAKER_MODULE_UTILS_TIME_H_
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#define NOERROR 0
#define ERROR 1
#define TIME_EXIT_OK 0
#define TIME_EXIT_NOTOK 1
#define TIME_MJD_START_IN_JD (double)(2400000.5)
#define TIME_IGREG_DATE_TO_JD (int)(588829)
#define TIME_IGREG_JD_TO_DATE (int)(2299161)
#define TIME_SEC_PER_MIN (float)(60.)
#define TIME_SEC_PER_HOUR (float)(3600.)
#define TIME_SEC_PER_DAY (float)(86400.)
#define TIME_25_MINS_IN_HOUR (float)(0.416667)
#define TIME_15_MINS_IN_HOUR (float)(0.25)
#define TIME_15_MINUTES_IN_MJD (double)(0.00976561)
#define TIME_25_MINUTES_IN_MJD (double)(0.017361)
#define TIME_30_MINUTES_IN_MJD (double)(0.0195312)
#define TIME_60_MINUTES_IN_MJD (double)(0.041667)
#define TIME_8_HOURS_IN_DAYS (double)(0.333333)
#define TIME_EPOCH_JD (double)(2436205.0)
#define TIME_EPOCH_MJD (double)(36204.5)
#define TIME_BIG_MJD (double)(TIME_EPOCH_MJD * 100.0)
typedef float fGEO_TIME_MJDate;
typedef double TGEO_TIME_MJDate;
typedef float fGEO_TIME_JDate;
typedef double TGEO_TIME_JDate;
struct SHORT_CDS_TIME
{
    unsigned short day;
    unsigned int msec;
};
inline double MJD_to_JD(double MJD){ return MJD + TIME_MJD_START_IN_JD; }
inline double JD_to_MJD(double JD){ return JD - TIME_MJD_START_IN_JD; }
inline void Days_to_Time(double days, int* seconds, int* minutes, int* hours)
{
   double dbSeconds,dbHours,dbMinutes;
   dbSeconds = (double)(TIME_SEC_PER_DAY) * days;
   dbHours = floor(dbSeconds / (double)(TIME_SEC_PER_HOUR));
   dbSeconds = dbSeconds - (double)(TIME_SEC_PER_HOUR) * dbHours;
   dbMinutes = floor(dbSeconds / (double)(TIME_SEC_PER_MIN) );
   dbSeconds = dbSeconds - (double)(TIME_SEC_PER_MIN) * dbMinutes;
  *seconds = (int)dbSeconds;
  *minutes = (int)dbMinutes;
  *hours = (int)dbHours;
}
inline int local_Midday_Hour(float longitude)
{
  if (longitude >= -180.0 && longitude <= 180.0)
     return (float) (12. - (longitude / 360. * 24.));
  else
  {
    return (float) (TIME_EXIT_NOTOK);
  }
}
inline int JD_to_Date(double JD, int* year, int* month, int* dayInMonth, int* hours, int* minutes, int* seconds)
{
  long ja,jalpha,jb,jc,jd,je;
  int iH,iM,iS;
  int iIJD;
  int days_in_prev_months[12] = {0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334};
  int iDayLeapYear;
  double dblFJD;
  iIJD = (int)floor(JD + 0.5);
 if (iIJD >= TIME_IGREG_JD_TO_DATE)
  {
  jalpha=(long)(((float) (iIJD-1867216)-0.25)/36524.25);
  ja=iIJD+1+jalpha-(long) (0.25*jalpha);
 } else if (iIJD < 0)
  {
  ja = iIJD + 36525 * (1-iIJD/36525);
 }
  else
  ja=iIJD;
 jb=ja+1524;
 jc=(long)(6680.0+((float) (jb-2439870)-122.1)/365.25);
 jd=(long)(365.25*jc);
 je=(long)((jb-jd)/30.6001);
 *dayInMonth=jb-jd-(long) (30.6001*je);
 *month=je-1;
 if (*month > 12) *month -= 12;
 *year=jc-4715;
 if (*month > 2) --(*year);
 if (*year <= 0) --(*year);
 if (iIJD < 0) *year -= 100*(1-iIJD/36525);
  dblFJD = JD - (double)(iIJD) + 0.5;
  Days_to_Time(dblFJD,&iS, &iM,&iH);
  *hours = iH;
  *minutes = iM;
  *seconds = iS;
  return TIME_EXIT_OK;
}
inline int Date_to_JD(int year, int month, int dayInMonth, int hours, int minutes, int seconds, double* JD)
{
  int iJulianMonth;
  int iIntYear,iJa;
  iIntYear = year;
 if (month > 2) {
  iJulianMonth = month + 1;
 } else
  {
  --iIntYear;
  iJulianMonth = month + 13;
 }
 *JD = (floor(365.25*iIntYear) + floor(30.6001*iJulianMonth) + (double)(dayInMonth) + 1720995.);
 if (dayInMonth+31L*(month+12L*year) >= TIME_IGREG_DATE_TO_JD)
  {
  iJa=(int)(0.01*iIntYear);
  *JD += 2-iJa+ (int)(0.25*iJa);
 }
  *JD += (seconds + (double)(TIME_SEC_PER_MIN) * minutes
         + (double)(TIME_SEC_PER_HOUR)*(hours - 12) ) / (double)(TIME_SEC_PER_DAY);
  return TIME_EXIT_OK;
}
inline int Date_to_MJD(int year, int month, int dayInMonth, int hours, int minutes, int seconds, double* MJD)
{
  int iStatus;
  double dbJd;
  iStatus = Date_to_JD(year, month, dayInMonth, hours, minutes, seconds, &dbJd);
  *MJD = JD_to_MJD(dbJd);
  return TIME_EXIT_OK;
}
inline int JDay_to_MonthDay(int year, int JDay, int* month, int* dayOfTheMonth)
{
  int iStatus;
  int iTmpMonth,iLoop,iTmpDay;
  int ArDays[] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
  int iDayLeapYear;
  int iDaycount;
  iDayLeapYear = 0;
  iTmpMonth = 1;
  iDaycount = 0;
  iTmpDay = 0;
  if ((( year%4 == 0 && year%100 != 0) || year%400 == 0))
  {
    iDayLeapYear = 1;
    ArDays[1] = 29;
  }
  for (iLoop=0;iLoop<12;iLoop++)
  {
    iTmpDay = JDay - iDaycount;
    iDaycount += ArDays[iLoop];
    if (iDaycount >= JDay) break;
    iTmpMonth++;
  }
  *month = iTmpMonth;
  *dayOfTheMonth = iTmpDay;
  return TIME_EXIT_OK;
}
inline int MJD_to_Time(double MJD, int* yearOut, int* dayOfTheYear, int* hourOut, int* minOut)
{
  int iYear,iMonth,iDayofTheMounth;
  int iHours,iMinutes,iSeconds;
  int days_in_prev_months[12] = {0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334};
  int iDayLeapYear;
  double dbJd;
  int iStatus;
  iDayLeapYear = 0;
  dbJd = MJD_to_JD(MJD);
  iStatus = JD_to_Date(dbJd, &iYear, &iMonth, &iDayofTheMounth, &iHours, &iMinutes, &iSeconds);
  if ( (( iYear%4 == 0 && iYear%100 != 0) || iYear%400 == 0) && (iMonth > 2) ) iDayLeapYear = 1;
  *dayOfTheYear = days_in_prev_months[iMonth-1] + iDayLeapYear + iDayofTheMounth;
  *yearOut = iYear;
  *hourOut = iHours;
  *minOut = iMinutes;
  return TIME_EXIT_OK;
}
inline int JD_to_DateTime(double JD, int* year, int* dayInYear, int* month, int* dayInMonth, int* hours, int* minutes, int* seconds)
{
  long ja,jalpha,jb,jc,jd,je;
  int iH,iM,iS;
  int iIJD;
  int days_in_prev_months[12] = {0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334};
  int iDayLeapYear;
  double dblFJD;
  iIJD = (int) floor(JD + 0.5);
 if (iIJD >= TIME_IGREG_JD_TO_DATE)
  {
  jalpha=(long)(((float) (iIJD-1867216)-0.25)/36524.25);
  ja=iIJD+1+jalpha-(long) (0.25*jalpha);
 } else if (iIJD < 0)
  {
  ja = iIJD + 36525 * (1-iIJD/36525);
 }
  else
  ja=iIJD;
 jb=ja+1524;
 jc=(long)(6680.0+((float) (jb-2439870)-122.1)/365.25);
 jd=(long)(365.25*jc);
 je=(long)((jb-jd)/30.6001);
 *dayInMonth=jb-jd-(long) (30.6001*je);
 *month=je-1;
 if (*month > 12) *month -= 12;
 *year=jc-4715;
 if (*month > 2) --(*year);
 if (*year <= 0) --(*year);
 if (iIJD < 0) *year -= 100*(1-iIJD/36525);
  dblFJD = JD - (double)(iIJD) + 0.5;
  Days_to_Time(dblFJD, &iS, &iM, &iH);
  *hours = iH;
  *minutes = iM;
  *seconds = iS;
  iDayLeapYear = 0;
  if ((( *year%4 == 0 && *year%100 != 0) || *year%400 == 0) && (*month > 2) ) iDayLeapYear = 1;
  *dayInYear = days_in_prev_months[*month-1] + iDayLeapYear + *dayInMonth;
  return TIME_EXIT_OK;
}
inline int SHRCDSTIME_Compare(SHORT_CDS_TIME* time1, SHORT_CDS_TIME* time2)
{
    if (time1 -> day > time2 -> day)
        return 1;
    else if (time1 -> day < time2 -> day)
        return -1;
    else
    {
        if (time1 -> msec > time2 -> msec)
            return 1;
        else if (time1 -> msec < time2 -> msec)
            return -1;
        else
            return 0;
    }
    return 0;
}
inline int CDS_TIME_to_DATE(SHORT_CDS_TIME *pstCDSTime,
                                       int *iYear,
                                       int *iMonth,
                                       int *iDay,
                                       float *flHour,
                                       float *flMin,
                                       float *flSec)
{
    int iStatus;
    int iRet;
    double MJDJulianDay;
    int iCurrentYear,
        iCurrentMonth,
        iCurrentDay,
        iCurrentrHours,
        iCurrentMin,
        iCurrentSec;
    float flMilli,flRest;
    iStatus = NOERROR;
    MJDJulianDay = TIME_EPOCH_JD + pstCDSTime->day;
    iRet = JD_to_Date(MJDJulianDay, &iCurrentYear, &iCurrentMonth, &iCurrentDay, &iCurrentrHours, &iCurrentMin, &iCurrentSec);
    if (iRet != NOERROR) iStatus = iRet;
    *iYear = iCurrentYear;
    *iMonth = iCurrentMonth;
    *iDay = iCurrentDay;
    flMilli = pstCDSTime->msec;
    *flHour = (float)floor(flMilli/(3600000L));
    flRest = flMilli - (*flHour)*3600000;
    *flMin = (float)floor(flRest/60000L);
    flRest = flRest - (*flMin)*60000;
    *flSec = (float)floor(flRest/1000);
    return iStatus;
}
#endif

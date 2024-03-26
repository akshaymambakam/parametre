%option noyywrap

%{
#include <stdio.h>

#define YY_DECL int yylex()

#include "ptre.tab.h"

extern int max_dot;

%}

%%

[ \t]	; // ignore all whitespace
[0-9]+\.[0-9]+ 	{yylval.fval = atof(yytext);
					int i=0;
					int dotp;
					while(yytext[i]!=0){
						if(yytext[i]=='.'){
							dotp=i;
						}
						i++;
					}
					max_dot = (i-dotp)>max_dot?(i-dotp):max_dot;
					return T_FLOAT;}
[0-9]+		   {yylval.ival = atoi(yytext); return T_INT;}
x[0-9]+   	   {
                 yylval.ival = atoi(yytext+1);
                 return X_STRING;
               }
y[0-9]+   	   {
                 yylval.ival = atoi(yytext+1);
                 return Y_STRING;
               }
z[0-9]+        {
                 yylval.ival = atoi(yytext+1);
                 return Z_STRING;
               }
p[0-9]+        {
                 yylval.ival = atoi(yytext+1);
                 return P_STRING;
               }

[ \t\r\n]+     { ; }

"^+"    {return K_PLUS;}
"+"		{return T_PLUS;}
"-"		{return T_MINUS;}
"*"		{return T_MULTIPLY;}
"<="    {yylval.ival=1; return LEQ;}
"<"     {yylval.ival=2; return LE;}
">="    {yylval.ival=3; return GEQ;}
">"     {yylval.ival=4; return GE;}
"MAX"   {yylval.ival=1; return T_MAX;}
"MIN"   {yylval.ival=2; return T_MIN;}
"DIFF"  {yylval.ival=3; return T_DIFF;}
"/"		{return T_DIVIDE;}
"&"    {return T_AND;}
"|"    {return T_OR;}
"%"    {return T_SEQCOMP;}
"["    {return T_SLEFT;}
"]"    {return T_SRIGHT;}
"("		{return T_LEFT;}
")"		{return T_RIGHT;}
"!"     {return T_NOT;}
"?"     {return T_PRE;}
"$"     {return T_POST;}
"EPS"   {return T_EPS;}
"+@"    {yylval.ival = 1; return T_RISE;}
"-@"     {yylval.ival = 0; return T_FALL;}

%%
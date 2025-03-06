%%
%% This is file `RUCPresentation.cls',
%% generated with the docstrip utility.
%%
%% The original source files were:
%%
%% classesOfRUC.dtx  (with options: `basis,pres')
%% classesOfRUC.dtx  (with options: `end')
%% ----------------------------------------------------------------
%% classesOfRUC ---
%% 
%% E-mail: anakin@ruc.dk
%% Released under the LaTeX Project Public License v1.3c or later
%% See http://www.latex-proj.org/lppl.txt
%% ----------------------------------------------------------------
\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{RUCPresentation}
   [2024/01/20 v1.0 Collection of RUC LaTeX classes by Anakin, anakin@ruc.dk]
%% ----------------------------------------------------------------
%% Utility Packages
%% ----------------------------------------------------------------
\RequirePackage{etoolbox}
\RequirePackage{iftex}
\RequirePackage{xifthen}
\RequirePackage{xkeyval}
%% ----------------------------------------------------------------
%% Class Tags
%% ----------------------------------------------------------------
\let\tag@RUC@Presentation@\iftrue\relax
\newif\iftagDoubleCol@
\newif\iftagSource@
\newif\iftagBibble@
\newif\iftagGlossel@
\newif\iftagFinal@
%% ----------------------------------------------------------------
%% Class Options
%% ----------------------------------------------------------------
\DeclareOption{10pt}{\PassOptionsToClass{\CurrentOption}{report}}
\DeclareOption{11pt}{\PassOptionsToClass{\CurrentOption}{report}}
\DeclareOption{12pt}{\PassOptionsToClass{\CurrentOption}{report}}
\DeclareOption{oneside}{
    \PassOptionsToClass{\CurrentOption}{report}
    \PassOptionsToPackage{centering}{geometry}
}
\DeclareOption{twoside}{\PassOptionsToClass{\CurrentOption}{report}}
\DeclareOption{a4paper}{
    \PassOptionsToClass{\CurrentOption}{report}
    \PassOptionsToPackage{\CurrentOption}{geometry}
}
\DeclareOption{aspectratio=169}{\PassOptionsToClass{\CurrentOption}{beamer}}
\newif\if@neverindent\@neverindentfalse
\DeclareOption{noindent}{\@neverindenttrue}
\DeclareOption{sOuRcEs}{\tagSource@true}
\DeclareOption{printBibble}{\tagBibble@true}
\DeclareOption{printGloss}{\tagGlossel@true}
\ExecuteOptions{11pt}
\ProcessOptions\relax
\LoadClass{beamer}
%% ----------------------------------------------------------------
%% Font Specifications
%% ----------------------------------------------------------------
\ifTUTeX
    % \usepackage[utf8]{luainputenc}
    % \usepackage[T1]{fontenc}
    % \usepackage{pxfonts}
    \RequirePackage[no-math]{fontspec}
    \RequirePackage[MnSymbol]{mathspec}
    \RequirePackage{unicode-math}
    % \defaultfontfeatures{Ligatures={TeX}}
    \setmainfont{Minion Pro}
\fi
\ifPDFTeX
    \RequirePackage[utf8]{inputenc} % utf8-unicode compatible input
    \RequirePackage[T1]{fontenc} % Font type specifier
\fi
\ifPDFTeX
    \RequirePackage[draft]{microtype}
\fi
\ifPDFTeX
    \PassOptionsToPackage{pdftex}{graphicx}
    \PassOptionsToPackage{pdftex}{hyperref}
\fi
\RequirePackage{amsmath,mathtools}
\RequirePackage{cancel}  % Striking through text \cancel{x}
\RequirePackage{geometry}
\RequirePackage[table, dvipsnames]{xcolor}
\RequirePackage{graphicx}
\RequirePackage{float}
\RequirePackage{caption,subcaption}
%% TODO Figure out if these fuck up beamer
% \RequirePackage{enumitem}
% \RequirePackage{fancyhdr} % Allows for easy fancy page layouts
% \RequirePackage{titlesec} % Allows for easy editing of of headers
\RequirePackage{listings}
\RequirePackage{setspace} % Line spacing
%% ----------------------------------------------------------------
%% TeX primitive redefinitions
%% ----------------------------------------------------------------
\tolerance=9999
\emergencystretch=10\p@
\hyphenpenalty=10000
\exhyphenpenalty=100
\interlinepenalty=10000
\let\@afterindentfalse\@afterindenttrue % Stolen from noindent package
\@afterindenttrue
\ifx\if@neverindent\iftrue
\setlength\parindent\z@
\else
    \setlength\parindent{1em}
\fi
\setlength\parskip{5\p@ plus\tw@\p@ minus\tw@\p@}
\setlength\marginparwidth{368\p@}
\setlength\headheight{14\p@}
%% ----------------------------------------------------------------
%% Microtype
%% ----------------------------------------------------------------
\ifPDFTeX
    \microtypesetup{%
        activate={true,nocompatibility}, %
        tracking=true, kerning=true,%
        spacing=true, factor=1100,%
        stretch=10, shrink=10}
    \microtypecontext{spacing=nonfrench}
    \SetProtrusion{encoding={*}, family={*}, series={*}, size={6,7}}{%
        1={ ,750}, 2={ ,500}, 3={ ,500}, 4={ ,500}, 5={ ,500},%
        6={ ,500}, 7={ ,600}, 8={ ,500}, 9={ ,500}, 0={ ,500}}
    \SetExtraKerning[unit=space]{%
        encoding={*}, family={*}, series={*}, %
        size={scriptsize,footnotesize,small,normalsize}}{%
        \textendash={400,400}, % en-dash, add more space around it
        "28={ ,150}, % left bracket, add space from right
        "29={150, }, % right bracket, add space from left
        \textquotedblleft={ ,150}, %
        \textquotedblright={150, }, %
        \textquoteleft={ ,150}, %
        \textquoteright={150, }}%
    \SetExtraKerning[unit=space] % For titles
       {encoding={*}, family={qhv}, series={b}, size={large,Large}}
       {1={-200,-200}, \textendash={400,400}}
    \SetTracking{encoding={*}, shape=sc}{40} % Fucks with small caps
\fi
%% ----------------------------------------------------------------
%% Colors
%% ----------------------------------------------------------------
\definecolor{myred}{rgb}{0.8,0.0,0.0}
\definecolor{myblue}{rgb}{0.8,0.0,0.0}
\definecolor{codebackgrounds}{rgb}{0.95,0.95,0.95}
\definecolor{codecomments}{rgb}{0.6,0.6,0.6}
\definecolor{codeidentifiers}{rgb}{0.4,0.0,0.4}
\definecolor{codekeywords}{rgb}{0.7,0.5,0.1}
\definecolor{codenumbers}{rgb}{0.5,0.5,0.5}
\definecolor{codestrings}{rgb}{0.22,0.67,0.22}
%% ----------------------------------------------------------------
%% Figure and Table formatting
%% ----------------------------------------------------------------
%% ----------------------------------------------------------------
%% Class Option Specifications
%% ----------------------------------------------------------------
\iftagSource@
    \PassOptionsToPackage{%
        sorting     = none,%
        maxbibnames = 99,%
        url         = false,%
        style       = numeric-comp,%
    }{biblatex}
    \PassOptionsToPackage{%
        nonumberlist,%
        toc, %
        xindy, %
        style   = index, %
    }{glossaries}
    \AtEndOfClass{
        \glsdisablehyper
        \setacronymstyle{long-sc-short}
    }
\fi
%% ----------------------------------------------------------------
%% Powerful Tools
%% ----------------------------------------------------------------
% Defining a new coordinate system for the page:
%
% --------------------------
% |(-1,1)    (0,1)    (1,1)|
% |                        |
% |(-1,0)    (0,0)    (1,0)|
% |                        |
% |(-1,-1)   (0,-1)  (1,-1)|
% --------------------------
\AtEndPreamble{
    \def\parsecomma#1,#2\endparsecomma{\def\page@x{#1}\def\page@y{#2}}
    \tikzdeclarecoordinatesystem{frame}{
        \parsecomma#1\endparsecomma
        \pgfpointanchor{current page}{north east}
        % Save the upper right corner
        \pgf@xc=\pgf@x%
        \pgf@yc=\pgf@y%
        % save the lower left corner
        \pgfpointanchor{current page}{south west}
        \pgf@xb=\pgf@x%
        \pgf@yb=\pgf@y%
        % Transform to the correct placement
        \pgfmathparse{(\pgf@xc-\pgf@xb)/2.*\page@x+(\pgf@xc+\pgf@xb)/2.}
        \expandafter\pgf@x\expandafter=\pgfmathresult pt
        \pgfmathparse{(\pgf@yc-\pgf@yb)/2.*\page@y+(\pgf@yc+\pgf@yb)/2.}
        \expandafter\pgf@y\expandafter=\pgfmathresult pt
    }
}
%
% \begin{tikzpicture}[remember picture, overlay]
%     \node (image) at (frame cs:0.5,0) {
%             \includegraphics[width=0.5\framewidth]
%             {image/path.file}
%         };
% \end{tikzpicture}
%
%
%
%% ----------------------------------------------------------------
%% House Keeping
%% ----------------------------------------------------------------
\author{Anakin}
\endinput
%% File: classesOfRUC.dtx
%% 
%% Copyright (C) 2024 by Anakin anakin@ruc.dk
%% -----------------------------------------------------------
%% 
%% This work may be distributed and/or modified under the
%% conditions of the LaTeX Project Public License (LPPL), either
%% version 1.3c of this license or (at your option) any later
%% version. The latest version of this license is in the file:
%% 
%%     http://www.latex-proj.org/lppl.txt
%% 
%% This work is "maintained" (as per LPPL maintenance status) by
%% Anakin.
%% 
%% This work consists of the file  classesOfRUC.dtx
%% and the derived files           classesOfRUC.ins,
%%                                 classesOfRUC.pdf,
%%                                 RUCAssignment.cls,
%%                                 RUCDSA.cls,
%%                                 RUCFSIM.cls,
%%                                 RUCSTAT.cls,
%%                                 RUCAECB.cls,
%%                                 RUCBIOINF.cls,
%%                                 RUCProject.cls, and
%%                                 RUCPresentation.cls
%% 
%%
%% End of file `RUCPresentation.cls'.

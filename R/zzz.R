.onLoad <- function(libname, pkgname) {

  #Import BED file
  data("bedfileMito")

  check_packages()
  checkDownloads()
  message("Installation DONE")

  checkHMTVAR_online()
  return(NA)
}

check_packages <- function() {
  if (!requireNamespace("ExomeDepth", quietly = TRUE)) {
    remotes::install_cran("ExomeDepth")
    install.packages("ExomeDepth", dependencies = TRUE)
  }
  library(ExomeDepth)
  if (!requireNamespace("openxlsx", quietly = TRUE)) {
    remotes::install_cran("openxlsx")
  }
  library(openxlsx)
  if (!requireNamespace("stringr", quietly = TRUE)) {
    remotes::install_cran("stringr")
  }
  library(stringr)
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    remotes::install_cran("ggplot2")
  }
  library(ggplot2)
  if (!requireNamespace("rjson", quietly = TRUE)) {
    remotes::install_cran("rjson")
  }
  library(rjson)
  if (!requireNamespace("rstudioapi", quietly = TRUE)) {
    remotes::install_cran("rstudioapi")
  }
  library(rstudioapi)
  if (!requireNamespace("tibble", quietly = TRUE)) {
    remotes::install_cran("tibble")
  }
  library(tibble)
  if (!requireNamespace("tidyverse", quietly = TRUE)) {
    remotes::install_cran("tidyverse")
  }
  library(tidyverse)
  if (!requireNamespace("showtext", quietly = TRUE)) {
    remotes::install_cran("showtext")
  }
  library(showtext)
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    remotes::install_cran("dplyr")
  }
  library(dplyr)
  if (!requireNamespace("ggtext", quietly = TRUE)) {
    remotes::install_cran("ggtext")
  }
  library(ggtext)
  if (!requireNamespace("readxl", quietly = TRUE)) {
    remotes::install_cran("readxl")
  }
  library(readxl)
  if (!requireNamespace("magrittr", quietly = TRUE)) {
    remotes::install_cran("magrittr")
  }
  library(magrittr)
  if (!requireNamespace("httr", quietly = TRUE)) {
    remotes::install_cran("httr")
  }
  library(httr)
  if (!requireNamespace("Rsamtools", quietly = TRUE)) {
    remotes::install_cran("Rsamtools")
  }
  library(Rsamtools)
  message("The R packages required have been successfully installed")
}

checkRequirements <- function() {
  OS <- tolower(system2("lsb_release", "-d", stdout = TRUE))
  # For Ubuntu
  if (grepl("ubuntu", OS) || grepl("kali", OS) ||
      grepl("mint", OS) || grepl("oracle", OS)){
    needed_packages_ubuntu <- c("bzip2", "libncurses5-dev", "unzip", "rpm2cpio", "cpio", "gzip", "tar", "wget", "java")
    for (i in needed_packages_ubuntu) {
      reply <- paste(system2("dpkg", sprintf("-s %s", i), stdout = TRUE), collapse = " ")
      if (grepl("Status: install ok installed", reply, ignore.case = TRUE)) {
        index <- grep(sprintf("%s", i), needed_packages_ubuntu)
        needed_packages_ubuntu[-index]
      }
    }
    # For RedHat
  } else if (grepl("redhat", OS) || grepl("fedora", OS) ||
             grepl("centos", OS)) {
    needed_packages_rpm <- c("bzip2", "ncurses-devel", "unzip", "rpm2cpio", "cpio", "gzip", "tar", "wget", "java")
    for (i in needed_packages_rpm) {
      reply <- paste(system2("rpm", sprintf("-q %s", i), stdout = TRUE), collapse = " ")

      if (!grepl("is not installed", reply, ignore.case = TRUE)) {
        index <- grep(sprintf("%s", i), needed_packages_rpm)
        needed_packages_rpm[-index]
      }
      needed_packages <- needed_packages_rpm
    }
  }
  #return(needed_packages)
  return(NULL)
}


# BWA
#' @title downloadBWA
#' @description Downloads and decompresses the BWA software
#' @return The path where the .exe file is
downloadBWA <- function(mitor_sof) {
  dir.create(sprintf("%s/BWA", mitor_sof))

  bwa_url1 <- "https://download.opensuse.org/repositories/home:/vojtaeus/15.4/x86_64/bwa-0.7.17-lp154.6.1.x86_64.rpm"
  bwa_dir1 <- file.path(mitor_sof, "BWA")
  #system2("wget", bwa_rpm_url, wait = TRUE, stdout = NULL, stderr = NULL, cwd = bwa_dir1)
  system2("wget", args = c(bwa_url1, "-P", bwa_dir1), wait = TRUE, stdout = NULL, stderr = NULL)
  #system2("rpm2cpio", "bwa-0.7.17-lp154.6.1.x86_64.rpm | cpio -id", wait = TRUE)
  system2("rpm2cpio", sprintf("%s/bwa-0.7.17-lp154.6.1.x86_64.rpm | cpio -D %s -idmv", bwa_dir1, bwa_dir1), wait = TRUE)


  dir.create(sprintf("%s/BWA/bwa-0.7.17-lp154.6.1.src", mitor_sof))
  #setwd(sprintf("%s/BWA/bwa-0.7.17-lp154.6.1.src", mitor_sof))
  bwa_url2 <-"https://download.opensuse.org/repositories/home:/vojtaeus/15.4/src/bwa-0.7.17-lp154.6.1.src.rpm"
  bwa_dir2 <- file.path(mitor_sof, "BWA/bwa-0.7.17-lp154.6.1.src")
  #system2("wget", bwa_url2, wait = TRUE, stdout = NULL, stderr = NULL)
  system2("wget", args = c(bwa_url2, "-P", bwa_dir2), wait = TRUE, stdout = NULL, stderr = NULL)
  #system2("rpm2cpio", "bwa-0.7.17-lp154.6.1.src.rpm | cpio -i", wait=TRUE)
  system2("rpm2cpio", sprintf("%s/bwa-0.7.17-lp154.6.1.src.rpm | cpio -D %s -idmv", bwa_dir2, bwa_dir2), wait = TRUE)


  #setwd(sprintf("%s/BWA", mitor_sof))
  bwa_url3 <-"https://download.opensuse.org/repositories/home:/vojtaeus/15.4/i586/bwa-0.7.17-lp154.6.1.i586.rpm"
  #system2("wget", bwa_url3, wait = TRUE, stdout = NULL, stderr = NULL)
  system2("wget", args = c(bwa_url3, "-P", bwa_dir1), wait = TRUE, stdout = NULL, stderr = NULL)
  #system2("rpm2cpio", "bwa-0.7.17-lp154.6.1.i586.rpm | cpio -id", wait=TRUE)
  system2("rpm2cpio", sprintf("%s/bwa-0.7.17-lp154.6.1.i586.rpm | cpio -D %s -idmv", bwa_dir1, bwa_dir1), wait = TRUE)

  return(sprintf('%s/BWA/usr/bin/bwa', mitor_sof))
}

# GATK
#' @title downloadGATK
#' @description Downloads and decompresses the GATK software
#' @return The path where the .exe file is located
downloadGATK <- function(mitor_sof) {
  #setwd(mitor_sof)
  URL <- 'https://github.com/broadinstitute/gatk/releases/download/4.3.0.0/gatk-4.3.0.0.zip'
  #system2("wget", URL, wait = TRUE)
  system2("wget", args = c(URL, "-P", mitor_sof), wait = TRUE, stdout = NULL, stderr = NULL)
  file <- basename(URL)
  file_dir <- file.path(mitor_sof, file)
  filedc <- substr(file, start = 0, stop = (nchar(file) - 4))
  #system2("unzip", sprintf(file), wait = TRUE)
  system2(sprintf("unzip %s -d %s", file_dir, mitor_sof), wait = TRUE)
  unzip(zipfile = file_dir, exdir = mitor_sof)
  return(sprintf('%s/%s/gatk-package-4.3.0.0-local.jar', mitor_sof, filedc))
}

# PICARD
#' @title downloadPICARD
#' @description Downloads and decompresses the PICARD software
#' @return The path where the .exe file is located
downloadPICARD <- function(mitor_sof) {
  URL <- "https://github.com/broadinstitute/picard/archive/refs/tags/2.27.5.tar.gz"
  #system2("wget", URL)
  system2("wget", args = c(URL, "-P", mitor_sof), wait = TRUE, stdout = NULL, stderr = NULL)
  system2("gzip" , sprintf("-d %s/2.27.5.tar.gz", mitor_sof))
  system2("tar" , sprintf("-xvf %s/2.27.5.tar -C %s", mitor_sof, mitor_sof))
  #system2("rm", substr(basename(URL), start = 0, stop = (nchar(basename(URL)) - 3)))
  file.remove(sprintf("%s/2.27.5.tar", mitor_sof))

  #setwd(sprintf("%s/picard-2.27.5", mitor_sof))
  #system2("wget", "https://github.com/broadinstitute/picard/releases/download/2.27.5/picard.jar")
  URL2 <- "https://github.com/broadinstitute/picard/releases/download/2.27.5/picard.jar"
  dir2 <- sprintf("%s/picard-2.27.5", mitor_sof)
  system2("wget", args = c(URL2, "-P", dir2), wait = TRUE, stdout = NULL, stderr = NULL)

  return(sprintf('%s/picard-2.27.5/picard.jar', mitor_sof))
}

# SAMTOOLS
#' @title downloadSamtools
#' @description Downloads and decompresses the Samtools software
#' @return The path where the .exe file is located
downloadSamtools <- function(mitor_sof) {
  dir.create(sprintf("%s/Samtools", mitor_sof))
  dir <- sprintf("%s/Samtools", mitor_sof)
  URL <- "https://github.com/samtools/samtools/releases/download/1.16.1/samtools-1.16.1.tar.bz2"
  system2("wget", args = c(URL, "-P", dir), wait = TRUE, stdout = NULL, stderr = NULL)

  #system2("bzip2", c("-d", basename(URL)), wait = TRUE)
  system2("bzip2", sprintf("-d %s/Samtools/samtools-1.16.1.tar.bz2", mitor_sof))
  system2("tar", c("-xvf", sprintf("%s/Samtools/%s -C %s", mitor_sof, list.files(sprintf("%s/Samtools", mitor_sof)), dir)))

  #file.remove(list.files(sprintf("%s/Samtools", mitor_sof))[stringr::str_detect(list.files(sprintf("%s/Samtools", mitor_sof)), ".tar")])
  file.remove(sprintf("%s/samtools-1.16.1.tar", dir))

  #setwd(sprintf("%s", list.files(sprintf("%s/Samtools", mitor_sof))))
  #system2("./configure")
  #system2("make")

  samtools_dir <- sprintf("%s/samtools-1.16.1", dir)
  system(paste("cd", shQuote(samtools_dir), "&& ./configure"))
  system(paste("cd", shQuote(samtools_dir), "&& make"))

  return(sprintf("%s/Samtools/samtools-1.16.1/samtools", mitor_sof))
}


#' @title downloadHG38
#' @description Downloads and decompresses the HG38 mitochodnrial DNA reference.
#' @return The path where the .fasta reference file is located
downloadHG38 <- function(mitor_sof) {
  dir.create(sprintf("%s/RefHG38", mitor_sof))
  URL <- "https://ftp.ensembl.org/pub/release-108/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.MT.fa.gz"
  dir <- sprintf("%s/RefHG38", mitor_sof)
  system2("wget", args = c(URL, "-P", dir), wait = TRUE, stdout = NULL, stderr = NULL)
  system2("gzip" , sprintf("-d %s/%s", dir, basename(URL)))

  reference <- substr(basename(URL), 1, nchar(basename(URL)) - 3)
  referenceN <- stringr::str_replace(reference, "fa", "fasta")
  referenceN <- sprintf("%s/%s", dir, referenceN)
  reference <- sprintf("%s/%s", dir, reference)

  system2("mv", sprintf("%s %s", reference, referenceN))

  BWA <- sprintf("%s/BWA/usr/bin/bwa", mitor_sof)
  GATK <- sprintf("%s/gatk-4.3.0.0/gatk-package-4.3.0.0-local.jar", mitor_sof)
  SAMTOOLS <- sprintf("%s/Samtools/samtools-1.16.1/samtools", mitor_sof)

  system2(BWA, sprintf("index -a is %s", referenceN))
  system2("java", sprintf("-jar %s CreateSequenceDictionary -R %s", GATK, referenceN))
  system2(SAMTOOLS, sprintf("faidx %s", referenceN))
  return(referenceN)
}

#' @title downloadMitoRSoftwares
#' @description Downloads all the necessary softwares and the HG38 file reference to use the MitoR software.
#' @return After each software is correctly downloaded and installed, a message from each one will be displayed on the terminal.
#' In case it occurs an error because of a missing package on your Linux OS, a message will be displayed on the terminal letting you know which packages must be available to download the software.
downloadMitoRSoftwares <- function(mitor_sof) {

  libPath <- dirname(system.file(package = "MitoR"))
  mitor_sof <- sprintf("%s/mitorDB/Softwares", libPath)

  required_packages <- checkRequirements()
  #required_packages <-
  if (!is.null(required_packages)) {
    stop(sprintf("The neccessary packages are not downloaded.
          Please install the following packages:
         %s

         On your terminal:

         for Ubuntu/Kali/Mint/Oracle:
            sudo apt-get install <package>

         for RedHat/CentOS:
            sudo -S yum install <package>

         for Fedora:
            sudo -S dnf install <package>

         for SUSE:
            sudo -S zypper install <package>


        Once you are done with the installation, please try loading MitoR again.
        ", paste(required_packages, collapse = ", ")))
    }
    tryCatch(
      expr = {
        SAMTOOLS <<- downloadSamtools(mitor_sof)
        system2(SAMTOOLS, "help")
      },
      error = function(e) {
        message("An error occured while performing the Samtools download and installation.
      Please remember that some packages are required and you must download them by yourself.
On the Linux command-line print:
------------------------------------------------------
    sudo apt install bzip2 libncurses5-dev tar
------------------------------------------------------
  libncurses5-dev for Debian or Ubuntu Linux or ncurses-devel for RPM-based Linux distributions")

        print(e)
      },

      finally = {
        message("-.Message from Samtools")
      }
    )
    tryCatch(
      expr = {
        GATK <<- downloadGATK(mitor_sof)
        system2("java", c(sprintf("-jar %s", GATK), "-h"))
        system2("rm", "~/gatk-4.3.0.0.zip")
      },
      error = function(e) {
        message("An error occured while performing the GATK download and installation.
      Please remember that some packages are required and you must download them by yourself.
On the Linux command-line print:
------------------------------------------------------
    sudo apt install unzip
------------------------------------------------------")
        print(e)
      },

      finally = {
        message("-.Message from GATK")
      }
    )
    tryCatch(
      expr = {
        BWA <<- downloadBWA(mitor_sof)
        system2(BWA)
        system2("rm", "~/bwa-0.7.17.tar")
      },
      error = function(e) {
        message("An error occured while performing the BWA download and installation.
      Please remember that some packages are required and you must download them by yourself.
On the Linux command-line print:
------------------------------------------------------
    sudo apt install wget rpm2cpio cpio
------------------------------------------------------")
        print(e)
      },
      warning = function(w) {
        message("An error occured while performing the BWA download and installation.
      Please remember that some packages are required and you must download them by yourself.
On the Linux command-line print:
------------------------------------------------------
    sudo apt install zip gzip rpm2cpio cpio
------------------------------------------------------")
        print(w)
      },
      finally = {
        message("-.Message from BWA")
      }
    )
    tryCatch(
      expr = {
        PICARD <<- downloadPICARD(mitor_sof)
        system2("java", sprintf("-jar %s -h", PICARD))
      },
      error = function(e) {
        message("An error occured while performing the PICARD download and installation.
      Please remember that some packages are required and you must download them by yourself.
On the Linux command-line print:
------------------------------------------------------
            sudo apt install gzip tar wget
------------------------------------------------------")
        print(e)
      },

      finally = {
        message("-.Message from PICARD")
      }
    )
    tryCatch(
      expr = {
        referencia <<- downloadHG38(mitor_sof)
        system2((sprintf("%s", BWA)), sprintf("index -a is %s", referencia))
      },
      error = function(e) {
        message("An error occured while performing the HG38 download and installation.
      Please remember that some packages are required and you must download them by yourself.
On the Linux command-line print:
------------------------------------------------------
            sudo apt install wget gzip java
------------------------------------------------------")
        print(e)
      },

      finally = {
        message("-.Message from NCBI")
      }
    )
  }

#' @title checkDownloads
#' @description Every time the mitor package is load, this function will be called just to make sure you have all the needed softwares to make it run.
#' In case there is one o more of them missing, they will automatically be downloaded.
#' @return When the download and installation is completed, a message will be displayed on the terminal

checkDownloads <- function() {

  libPath <- dirname(dirname(dirname(system.file(package = "MitoR"))))
  print(libPath)
  if (!(grepl("linux", libPath, ignore.case = TRUE))) {
    return(NULL)
  }

  if (!(file.exists(sprintf("%s/mitorDB/Softwares/BWA", libPath))) |
      !(file.exists(sprintf("%s/mitorDB/Softwares/Samtools", libPath))) |
      !(file.exists(sprintf("%s/mitorDB/Softwares/RefHG38", libPath)))
    ){
    message("The softwares will be downloaded. This process can take a few minutes.")
    dir.create(sprintf("%s/mitorDB", libPath))
    dir.create(sprintf("%s/mitorDB/Softwares", libPath))
    dir.create(sprintf("%s/mitorDB/DB", libPath))

    mitor_sof <- sprintf("%s/mitorDB/Softwares", libPath)
    mitor_db <- sprintf("%s/mitorDB/DB", libPath)

    #setwd(mitor_sof)
    downloadMitoRSoftwares(mitor_sof)

  } else {
    message("The softwares are already downloaded. We will make sure they can be used correctly.")
    mitor_sof <- sprintf("%s/mitorDB/Softwares", libPath)
    mitor_db <- sprintf("%s/mitorDB/DB", libPath)
    #setwd(mitor_sof)

    tryCatch(expr = {
      SAMTOOLS <- sprintf("%s/Samtools/samtools-1.16.1/samtools", mitor_sof)
      BWA <- sprintf("%s/BWA/usr/bin/bwa", mitor_sof)
      PICARD <- sprintf("%s/picard-2.27.5/picard.jar", mitor_sof)
      GATK <- sprintf("%s/gatk-4.3.0.0/gatk-package-4.3.0.0-local.jar", mitor_sof)
      referencia <- sprintf("%s/RefHG38/Homo_sapiens.GRCh38.dna.chromosome.MT.fasta", mitor_sof)
      system2(sprintf("%s", SAMTOOLS))
      system2(sprintf("%s", BWA))
      system2("java", sprintf("-jar %s -h", GATK))
      system2("java", sprintf("-jar %s -h", PICARD))
    },
    #warning = function(cond) {
    #  system2("rm", mitor_sof)
    #  dir.create(mitor_sof)
      #setwd(mitor_sof)
    #  downloadMitoRSoftwares(mitor_sof)
    #  message("We had to download everything again due to a missing file")
    #  return(NA)
    #},

    error = function(cond) {
      system2("rm", mitor_sof)
      dir.create(mitor_sof)
      #setwd(mitor_sof)
      downloadMitoRSoftwares(mitor_sof)
      message("We had to download everything again due to a missing file")
      return(NA)
    },
    finally = {
      message(" -.Download and Installation of softwares is done")
    }
    )
  }
}



checkHMTVAR_online <- function() {

  if (!httr::http_error(httr::GET("https://www.google.com"))) {
    libPath <- dirname(system.file(package = "MitoR"))
    mitor_db <- sprintf("%s/mitorDB/DB", libPath)

    if ("RDS_DB.rds" %in% list.files(mitor_db)) {
      RDS_DB <- readRDS(sprintf("%s/RDS_DB.rds", mitor_db))
      if (any(stringr::str_detect(RDS_DB[[1]]$Freq_HMTVAR , "Offline"))){
        Freq_MitoR <- RDS_DB[[1]]
        # Select columns that contain any value "Offline"
        valuesOffline <- Freq_MitoR[Freq_MitoR$Freq_HMTVAR == "Offline", ]
        valuesToUpdate <- buscarHMTVAR(valuesOffline)
        # Agregamos las filas modificadas en el mismo lugar que estaban antes las offline
        Freq_MitoR[valuesOffline, ] <- valuesToUpdate
        # Volvemos a guardar en el RDS_DB el dataframe actualizado y lo guardamos como RDS
        RDS_DB$Freq_MitoR <- Freq_MitoR
        saveRDS(RDS_DB, sprintf("%s/RDS_DB.rds", mitor_db))
      }
    }
  }
}



.onLoad <- function(libname, pkgname) {

  #Import BED file
  data("bedfileMito")

  check_packages()
  checkDownloads()
  message("Softwares DONE")

  checkHMTVAR_online()
  return(NA)
}

check_packages <- function() {
  #if (!requireNamespace("ExomeDepth", quietly = TRUE)) {
  #  remotes::install_cran("ExomeDepth")
  #}
  if (!requireNamespace("openxlsx", quietly = TRUE)) {
    remotes::install_cran("openxlsx")
  }
  if (!requireNamespace("stringr", quietly = TRUE)) {
    remotes::install_cran("stringr")
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    remotes::install_cran("ggplot2")
  }
  if (!requireNamespace("rjson", quietly = TRUE)) {
    remotes::install_cran("rjson")
  }
  if (!requireNamespace("rstudioapi", quietly = TRUE)) {
    remotes::install_cran("rstudioapi")
  }
  if (!requireNamespace("tibble", quietly = TRUE)) {
    remotes::install_cran("tibble")
  }
  if (!requireNamespace("tidyverse", quietly = TRUE)) {
    remotes::install_cran("tidyverse")
  }
  if (!requireNamespace("showtext", quietly = TRUE)) {
    remotes::install_cran("showtext")
  }
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    remotes::install_cran("dplyr")
  }
  if (!requireNamespace("ggtext", quietly = TRUE)) {
    remotes::install_cran("ggtext")
  }
  if (!requireNamespace("readxl", quietly = TRUE)) {
    remotes::install_cran("readxl")
  }
  if (!requireNamespace("magrittr", quietly = TRUE)) {
    remotes::install_cran("magrittr")
  }
  if (!requireNamespace("httr", quietly = TRUE)) {
    remotes::install_cran("httr")
  }
  #if (!requireNamespace("Rsamtools", quietly = TRUE)) {
  #  remotes::install_cran("Rsamtools")
  #}
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
      # A chequear que realmente devuelva esta afirmacion
      if (!grepl("is not installed", reply, ignore.case = TRUE)) {
        index <- grep(sprintf("%s", i), needed_packages_rpm)
        needed_packages_rpm[-index]
      }
      needed_packages <- needed_packages_rpm
    }
  }
  return(needed_packages)
}


# BWA
#' @title downloadBWA
#' @description Downloads and decompresses the BWA software
#' @return The path where the .exe file is
downloadBWA <- function(mitor_sof) {
  #setwd(mitor_sof)
  dir.create(sprintf("%s/BWA", mitor_sof))

  setwd(sprintf("%s/BWA", mitor_sof))
  system2("wget", "https://download.opensuse.org/repositories/home:/vojtaeus/15.4/x86_64/bwa-0.7.17-lp154.6.1.x86_64.rpm", wait=TRUE)
  system2("rpm2cpio", "bwa-0.7.17-lp154.6.1.x86_64.rpm | cpio -id", wait=TRUE)

  dir.create(sprintf("%s/BWA/bwa-0.7.17-lp154.6.1.src", mitor_sof))
  setwd(sprintf("%s/BWA/bwa-0.7.17-lp154.6.1.src", mitor_sof))
  system2("wget", "https://download.opensuse.org/repositories/home:/vojtaeus/15.4/src/bwa-0.7.17-lp154.6.1.src.rpm", wait=TRUE)
  system2("rpm2cpio", "bwa-0.7.17-lp154.6.1.src.rpm | cpio -i", wait=TRUE)

  setwd(sprintf("%s/BWA", mitor_sof))
  system2("wget", "https://download.opensuse.org/repositories/home:/vojtaeus/15.4/i586/bwa-0.7.17-lp154.6.1.i586.rpm", wait=TRUE)
  system2("rpm2cpio", "bwa-0.7.17-lp154.6.1.i586.rpm | cpio -id", wait=TRUE)

  return(sprintf('%s/BWA/usr/bin/bwa', mitor_sof))
}

# GATK
#' @title downloadGATK
#' @description Downloads and decompresses the GATK software
#' @return The path where the .exe file is located
downloadGATK <- function(mitor_sof) {
  #setwd(mitor_sof)
  URL <- 'https://github.com/broadinstitute/gatk/releases/download/4.3.0.0/gatk-4.3.0.0.zip'
  system2("wget", URL, wait = TRUE)
  file <- basename(URL)
  filedc <- substr(file,start = 0, stop = (nchar(file) - 4))
  system2("unzip", sprintf(file), wait = TRUE)
  setwd(sprintf("%s/%s", mitor_sof, filedc))
  return(sprintf('%s/%s/gatk-package-4.3.0.0-local.jar', mitor_sof, filedc))
}

# PICARD
#' @title downloadPICARD
#' @description Downloads and decompresses the PICARD software
#' @return The path where the .exe file is located
downloadPICARD <- function(mitor_sof) {
  #setwd(mitor_sof)
  URL <- "https://github.com/broadinstitute/picard/archive/refs/tags/2.27.5.tar.gz"
  system2("wget", URL)
  system2("gzip" , sprintf("-d %s/2.27.5.tar.gz", mitor_sof))
  system2("tar" , sprintf("-xvf %s/2.27.5.tar", mitor_sof))
  system2("rm", substr(basename(URL), start = 0, stop = (nchar(basename(URL)) - 3)))
  setwd(sprintf("%s/picard-2.27.5", mitor_sof))
  system2("wget", "https://github.com/broadinstitute/picard/releases/download/2.27.5/picard.jar")
  return(sprintf('%s/picard-2.27.5/picard.jar', mitor_sof))
}

# SAMTOOLS
#' @title downloadSamtools
#' @description Downloads and decompresses the Samtools software
#' @return The path where the .exe file is located
downloadSamtools <- function(mitor_sof) {
  #setwd(mitor_sof)
  dir.create(sprintf("%s/Samtools", mitor_sof))
  setwd(sprintf("%s/Samtools", mitor_sof))
  URL <- "https://github.com/samtools/samtools/releases/download/1.16.1/samtools-1.16.1.tar.bz2"
  system2("wget", URL, wait = TRUE)
  system2("bzip2", c("-d", basename(URL)), wait = TRUE)
  system2("tar" , c("-xvf", sprintf("%s/Samtools/%s", mitor_sof, list.files(sprintf("%s/Samtools", mitor_sof)))))
  file.remove(list.files(sprintf("%s/Samtools", mitor_sof))[stringr::str_detect(list.files(sprintf("%s/Samtools", mitor_sof)), ".tar")])
  setwd(sprintf("%s", list.files(sprintf("%s/Samtools", mitor_sof))))
  system2("./configure")
  system2("make")
  return(sprintf("%s/Samtools/samtools-1.16.1/samtools", mitor_sof))
}

# REFERENCIA
#' @title downloadHG38
#' @description Downloads and decompresses the HG38 mitochodnrial DNA reference.
#' @return The path where the .fasta reference file is located
downloadHG38 <- function(mitor_sof) {
  #setwd(mitor_sof)
  dir.create(sprintf("%s/RefHG38", mitor_sof))
  setwd(sprintf("%s/RefHG38", mitor_sof))
  URL <- "https://ftp.ensembl.org/pub/release-108/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.MT.fa.gz"
  referencia <- substr(basename(URL), 1, nchar(basename(URL)) - 3)
  system2("wget", URL, wait = TRUE)
  system2("gzip" , sprintf("-d %s", basename(URL)))
  referenciaN <- stringr::str_replace(referencia, "fa", "fasta")
  system2("mv", sprintf("%s %s", referencia, referenciaN))
  #chmod +x /home/daniela/MitoRSoftware/BWA/usr/bin/bwa

  BWA <- sprintf("%s/BWA/usr/bin/bwa", mitor_sof)
  GATK <- sprintf("%s/gatk-4.3.0.0/gatk-package-4.3.0.0-local.jar", mitor_sof)
  SAMTOOLS <- sprintf("%s/Samtools/samtools-1.16.1/samtools", mitor_sof)

  system2(BWA, sprintf("index -a is %s", referenciaN))
  system2("java", sprintf("-jar %s CreateSequenceDictionary -R %s", GATK, referenciaN))
  system2(SAMTOOLS, sprintf("faidx %s", referenciaN))
  return(sprintf("%s/RefHG38/%s", mitor_sof, referenciaN))
}

#' @title downloadMitoRSoftwares
#' @description Downloads all the necessary softwares and the HG38 file reference to use the MitoR software.
#' @return After each software is correctly downloaded and installed, a message from each one will be displayed on the terminal.
#' In case it occurs an error because of a missing package on your Linux OS, a message will be displayed on the terminal letting you know which packages must be available to download the software.
downloadMitoRSoftwares <- function(mitor_sof) {

  libPath <- dirname(system.file(package = "MitoR"))
  mitor_sof <- sprintf("%s/mitorDB/Softwares", libPath)
  #setwd(mitor_sof)
  #wd <<- getwd()
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
      warning = function(w) {
        message("An error occured while performing the Samtools download and installation.
      Please remember that some packages are required and you must download them by yourself.
On the Linux command-line print:
------------------------------------------------------
    sudo apt install bzip2 libncurses5-dev tar
------------------------------------------------------
  libncurses5-dev for Debian or Ubuntu Linux distributions or ncurses-devel for RPM-based Linux distributions")
        print(w)
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
      warning = function(w) {
        message("An error occured while performing the GATK download and installation.
      Please remember that some packages are required and you must download them by yourself.
On the Linux command-line print:
------------------------------------------------------
    sudo apt install unzip
------------------------------------------------------")
        print(w)
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
      warning = function(w) {
        message("An error occured while performing the PICARD download and installation.
      Please remember that some packages are required and you must download them by yourself.
On the Linux command-line print:
------------------------------------------------------
            sudo apt install gzip tar wget
------------------------------------------------------")
        print(w)
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
      warning = function(w) {
        message("An error occured while performing the HG38 download and installation.
      Please remember that some packages are required and you must download them by yourself.
On the Linux command-line print:
------------------------------------------------------
            sudo apt install wget gzip java
------------------------------------------------------")
        print(w)
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

  if (!(file.exists(sprintf("%s/mitorDB", libPath)))) {
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



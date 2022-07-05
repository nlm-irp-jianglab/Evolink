install.packages('BiocManager',
                 dependencies='Depends',
                 repo = "http://cran.us.r-project.org")
# Roll back version to v0.1.8 would make loading ggtree work
install.packages("https://cran.r-project.org/src/contrib/Archive/rvcheck/rvcheck_0.1.8.tar.gz", repos = NULL)
BiocManager::install('ggtree')
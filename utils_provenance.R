add_git_info = function(filename){
  # TBD: make this function portable to operating systems other than Linux
  if(file.exists(".git")){
    git_repo = system("basename `git remote -v | head -1 | awk '{ print $(NF-1) }'` .git", intern = TRUE)
    
    git_info = system("git log --graph --all --decorate | head -3", intern = TRUE)
    
    commit = git_info[1]
    author = git_info[2]
    date = git_info[3]
    
    git_branch = system("git rev-parse --abbrev-ref HEAD", intern = TRUE)
    
    git_message = sprintf("#\n#\n### VERSION CONTROL INFORMATION ###\n# %s\n# %s\n# In Git repository %s, on branch %s,\n# %s", 
                          author, date, git_repo, git_branch, commit)
    
    write(git_message, filename, append = TRUE)
  }
}

add_R_session_info = function(filename){
  cat(c("#\n#\n### R SESSION INFORMATION: ###", 
        captureOutput(sessionInfo())), 
      sep = "\n# ", file = filename, append = TRUE)
}

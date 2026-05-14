

                                                                                                                                                  # 1. 重置到远程最新                                                                                                                              
git fetch origin                                                                                                                                 git reset --soft origin/master                                                                                                                                                                                                                                                                    # 2. 设置身份                                                                                                                                    
git config user.name "Xiaoting Hua"                                                                                                              git config user.email "xiaotinghua@zju.edu.cn"                                                                                                                                                                                                                                                    # 3. 添加所有新文件                                                                                                                              
git add -A                                    
                                           
     # 4. 以你的名字提交                                                                                                                              
git commit -m "feat: PWM-based pdif detection v2 + v1 CR pairing + DC validation + Hall degraded search"                                                                                                                                                                                          # 5. 强制推送（用 -f）                                                                                                                           
git push -f origin master    


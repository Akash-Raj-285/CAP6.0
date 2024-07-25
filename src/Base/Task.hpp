/* **********************************************************************
 * Copyright (C) 2019-2022, Claude Pruneau, Victor Gonzalez, Sumit Basu
 * All rights reserved.
 *
 * Based on the ROOT package and environment
 *
 * For the licensing terms see LICENSE.
 *
 * Author: Claude Pruneau,   04/01/2022
 *
 * *********************************************************************/
#ifndef CAP__Task
#define CAP__Task
#include "TClass.h"
#include "Aliases.hpp"
#include "Exceptions.hpp"
#include "ConfigurationManager.hpp"
#include "TaskAccountant.hpp"
#include "NamedObject.hpp"
#include "Manager.hpp"
#include "EnvironmentVariables.hpp"
#include "Timer.hpp"
#include "TFile.h"
#include <iostream>

namespace CAP
{

class Task
:
public ConfigurationManager,
public NamedObject,
public TaskAccountant,
public Timer,
public EnvironmentVariables
{
protected:
  //!
  //! Pointer to parent task if any
  //!
  Task * parent;

  //!
  //! Array of pointers to subTasks called by this task instance, once per event analyzed (or iteration generated by EventIterator task). If this instance carries out
  //! initialize, finalize, execute type operations, these are performed BEFORE the corresponding operations by the subTasks.
  //!
  std::vector<Task*> subTasks;

public:

  Task();
  Task(const Task & task);
  Task & operator=(const Task & task);
  virtual ~Task() {}

  virtual void setDefaultConfiguration();
  virtual void configure();
  virtual void initialize();
  virtual void execute();
  virtual void partial();
  virtual void finalize();
  virtual void reset();
  virtual void clear();
  virtual void print(std::ostream & output=std::cout,  
                     int style=0,
                     int size=50) const;

  void configureSubTasks();
  void initializeSubTasks();
  void executeSubTasks();
  void finalizeSubTasks();
  void partialSubTasks();
  void resetSubTasks();
  void clearSubTasks();
  void printSubTasks(std::ostream & output=std::cout,
                     int style=0,
                     int size=50) const;

  std::vector<Task*> & getSubTaks();
  const std::vector<Task*> & getSubTaks() const;
  bool hasSubTasks()  const;
  unsigned int getNSubTasks() const;
  Task * getSubTaskAt(unsigned int index);
  Task *  addSubTask(Task * task);

  String getParentTaskName() const;
  bool hasParent() const;
  Task * getParent() const;
  void setParent(Task * _parent);
  String getParentName() const;
  Task * getTaskAt(int depth);
  const Task * getTaskAt(int depth) const;
  String getReverseTaskPath(int depth=0) const;
  VectorString getTaskPathTokens() const;
  String getTaskPath(int depth=0) const;
  String getFullTaskPath() const;
  int getNAncestors() const;

  ClassDef(Task,0)
};

}

#endif /* CAP__Task */

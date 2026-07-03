# Software Design Document

## 1D-QM-Playground

**Version:** 0.1 (Draft)
**Date:** 2026-07-03
**Status:** Skeleton — sections pending

---

## Table of Contents

1. [Introduction](#1-introduction)
   - [1.1 Purpose of This Document](#11-purpose-of-this-document)
   - [1.2 Scope](#12-scope)
   - [1.3 Intended Audience](#13-intended-audience)
   - [1.4 Definitions, Acronyms, and Abbreviations](#14-definitions-acronyms-and-abbreviations)
   - [1.5 References](#15-references)
   - [1.6 Document Conventions](#16-document-conventions)
2. [System Overview](#2-system-overview)
   - [2.1 Background and Motivation](#21-background-and-motivation)
   - [2.2 Goals and Objectives](#22-goals-and-objectives)
   - [2.3 System Context](#23-system-context)
   - [2.4 Assumptions and Constraints](#24-assumptions-and-constraints)
3. [Requirements](#3-requirements)
   - [3.1 Functional Requirements](#31-functional-requirements)
   - [3.2 Non-Functional Requirements](#32-non-functional-requirements)
   - [3.3 Requirements Traceability Matrix](#33-requirements-traceability-matrix)
4. [Architectural Design](#4-architectural-design)
   - [4.1 Architectural Style and Rationale](#41-architectural-style-and-rationale)
   - [4.2 Component Overview](#42-component-overview)
   - [4.3 Data Flow](#43-data-flow)
   - [4.4 Execution View](#44-execution-view)
5. [Detailed Component Design](#5-detailed-component-design)
   - [5.1 Controller](#51-controller)
   - [5.2 TISE Solver](#52-tise-solver)
   - [5.3 TDSE Solver](#53-tdse-solver)
   - [5.4 Analysis Module](#54-analysis-module)
6. [Data Design](#6-data-design)
   - [6.1 Configuration Schema](#61-configuration-schema)
   - [6.2 Internal Data Structures](#62-internal-data-structures)
   - [6.3 Persistent Storage Format](#63-persistent-storage-format)
   - [6.4 Data Validation Rules](#64-data-validation-rules)
7. [Interface Design](#7-interface-design)
   - [7.1 External Interfaces](#71-external-interfaces)
   - [7.2 Inter-Component Interfaces](#72-inter-component-interfaces)
   - [7.3 API and Function Signatures](#73-api-and-function-signatures)
8. [Error Handling and Logging Strategy](#8-error-handling-and-logging-strategy)
9. [Testing Strategy](#9-testing-strategy)
   - [9.1 Unit Testing and TDD Approach](#91-unit-testing-and-tdd-approach)
   - [9.2 Integration Testing](#92-integration-testing)
   - [9.3 Verification and Validation](#93-verification-and-validation)
   - [9.4 Test Traceability](#94-test-traceability)
10. [Build, Configuration Management, and Deployment](#10-build-configuration-management-and-deployment)
    - [10.1 Build System and Dependencies](#101-build-system-and-dependencies)
    - [10.2 Directory Layout](#102-directory-layout)
    - [10.3 Version Control Conventions](#103-version-control-conventions)
11. [Appendices](#11-appendices)
    - [A. Glossary](#a-glossary)
    - [B. Open Design Questions](#b-open-design-questions)
    - [C. Revision History](#c-revision-history)

---

## 1. Introduction

### 1.1 Purpose of This Document

This Software Design Document (SDD) is the authoritative, single source of truth for the design of the 1D-QM-Playground system. Its purpose is to translate the physics goals and requirements of this project — a publication-quality B-spline based solver for the time-independent and time-dependent Schrödinger equations — into a concrete, unambiguous engineering plan that can be implemented, tested, and maintained with confidence. Members of this project bring deep expertise in physics and numerical methods, but may have less occasion to work with formal software engineering documents; this section accordingly explains why an SDD matters and how it will be used, so that its value is clear before the detailed technical content is filled in.

A design document exists first to separate the "what" from the "how" from the "when." Requirements — what the system must do — architecture — how the system is organized to do it — and implementation — the actual code — are three distinct concerns. Conflating them, by designing while coding, tends to produce systems where a change to one part has unpredictable effects on distant parts, because no one wrote down the boundaries between components. This document exists to make those boundaries explicit and durable before any of them are committed to code.

This separation is also what makes Test-Driven Development (TDD) possible. TDD is the practice of writing a test that specifies expected behavior before writing the code that implements it, which can only be done when the expected behavior is already known — that is, when the design's interfaces, inputs, outputs, and edge cases have been decided ahead of time. By fully specifying component interfaces and behaviors in this document, tests can be written first, so that "done" comes to mean "meets a written, testable specification" rather than "seems to work."

The same specificity enables interface-driven implementation. Once the contract between two pieces of the system — a function signature, a file format, a data schema — is fixed in this document, both sides of that contract can be implemented and tested independently, even in parallel, as long as each honors the agreed interface. For example, once the `config.yaml` schema is fixed here, the TISE solver, TDSE solver, and Analysis module can each be developed and unit-tested against that schema without waiting on each other or on the Controller.

This document is also the mechanism for requirements-driven design and implementation. Every design decision recorded here should trace back to a stated requirement, whether a physics capability, a performance need, or a usability need. This gives the project a way to check, at any point, whether a piece of code is necessary — does it trace to a requirement — and whether a requirement is satisfied — what code implements it. It equally gives a principled way to decline scope creep: a proposed feature that does not trace to a requirement does not belong in this version of the system.

Deciding these matters up front, and recording them where they can be reviewed, also reduces the accumulation of technical debt. Technical debt is the extra rework caused by choosing a quick, ad-hoc solution now instead of a well-considered one; the debt comes due later as bugs, as confusing code, or as work that must be redone under time pressure. Fixing a bad decision on paper, before code and tests depend on it, is far cheaper than fixing it afterward.

Taken together, these properties make the implementation process substantially more efficient. A contributor picking up any one component of this system — the TISE solver, say — should be able to read this document and know what the component is responsible for, what it receives as input and must produce as output, what other components depend on it, and how its correctness will be verified. This removes the need to reverse-engineer intent from code or to interrupt other contributors with clarifying questions, and it means the same document that guided implementation can later guide code review and onboarding.

This document will therefore be treated as living but authoritative: it should be updated whenever a design decision changes, and code should be expected to match what is written here. Where the two disagree, that disagreement is itself a defect to be resolved — either the document or the code is wrong, and both cases are worth fixing promptly.

### 1.2 Scope

### 1.3 Intended Audience

### 1.4 Definitions, Acronyms, and Abbreviations

### 1.5 References

### 1.6 Document Conventions

---

## 2. System Overview

### 2.1 Background and Motivation

### 2.2 Goals and Objectives

### 2.3 System Context

### 2.4 Assumptions and Constraints

---

## 3. Requirements

### 3.1 Functional Requirements

### 3.2 Non-Functional Requirements

### 3.3 Requirements Traceability Matrix

---

## 4. Architectural Design

### 4.1 Architectural Style and Rationale

### 4.2 Component Overview

### 4.3 Data Flow

### 4.4 Execution View

---

## 5. Detailed Component Design

### 5.1 Controller

#### 5.1.1 Responsibilities

#### 5.1.2 Interfaces

#### 5.1.3 Internal Design

#### 5.1.4 Error Handling

### 5.2 TISE Solver

#### 5.2.1 Responsibilities

#### 5.2.2 Interfaces

#### 5.2.3 Internal Design

#### 5.2.4 Error Handling

### 5.3 TDSE Solver

#### 5.3.1 Responsibilities

#### 5.3.2 Interfaces

#### 5.3.3 Internal Design

#### 5.3.4 Error Handling

### 5.4 Analysis Module

#### 5.4.1 Responsibilities

#### 5.4.2 Interfaces

#### 5.4.3 Internal Design

#### 5.4.4 Error Handling

---

## 6. Data Design

### 6.1 Configuration Schema

### 6.2 Internal Data Structures

### 6.3 Persistent Storage Format

### 6.4 Data Validation Rules

---

## 7. Interface Design

### 7.1 External Interfaces

### 7.2 Inter-Component Interfaces

### 7.3 API and Function Signatures

---

## 8. Error Handling and Logging Strategy

---

## 9. Testing Strategy

### 9.1 Unit Testing and TDD Approach

### 9.2 Integration Testing

### 9.3 Verification and Validation

### 9.4 Test Traceability

---

## 10. Build, Configuration Management, and Deployment

### 10.1 Build System and Dependencies

### 10.2 Directory Layout

### 10.3 Version Control Conventions

---

## 11. Appendices

### A. Glossary

### B. Open Design Questions

### C. Revision History

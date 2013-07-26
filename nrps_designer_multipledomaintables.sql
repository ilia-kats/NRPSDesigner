-- phpMyAdmin SQL Dump
-- version 3.2.4
-- http://www.phpmyadmin.net
--
-- Host: localhost
-- Erstellungszeit: 09. Juli 2013 um 14:35
-- Server Version: 5.1.44
-- PHP-Version: 5.3.1

SET SQL_MODE="NO_AUTO_VALUE_ON_ZERO";


/*!40101 SET @OLD_CHARACTER_SET_CLIENT=@@CHARACTER_SET_CLIENT */;
/*!40101 SET @OLD_CHARACTER_SET_RESULTS=@@CHARACTER_SET_RESULTS */;
/*!40101 SET @OLD_COLLATION_CONNECTION=@@COLLATION_CONNECTION */;
/*!40101 SET NAMES utf8 */;

--
-- Datenbank: `nrps_designer`
--
CREATE DATABASE `nrps_designer` DEFAULT CHARACTER SET utf8 COLLATE utf8_unicode_ci;
USE `nrps_designer`;

-- --------------------------------------------------------

--
-- Tabellenstruktur für Tabelle `domain_types_a`
--

CREATE TABLE IF NOT EXISTS `domain_types_a` (
  `pathway_id` int(20) NOT NULL,
  `domain_borders` varchar(20) COLLATE utf8_unicode_ci NOT NULL,
  `initiation_domain` tinyint(1) NOT NULL,
  `substrate_specificity_as_id` varchar(50) COLLATE utf8_unicode_ci NOT NULL,
  `uni_prot_id` varchar(20) COLLATE utf8_unicode_ci NOT NULL,
  `bb_id` varchar(20) COLLATE utf8_unicode_ci NOT NULL,
  `description` text COLLATE utf8_unicode_ci NOT NULL,
  `dna_sequence` text COLLATE utf8_unicode_ci NOT NULL,
  `native_linker_before` text COLLATE utf8_unicode_ci NOT NULL,
  `native_linker_after` text COLLATE utf8_unicode_ci NOT NULL,
  PRIMARY KEY (`pathway_id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8 COLLATE=utf8_unicode_ci;

--
-- Daten für Tabelle `domain_types_a`
--


-- --------------------------------------------------------

--
-- Tabellenstruktur für Tabelle `domain_types_a_c`
--

CREATE TABLE IF NOT EXISTS `domain_types_a_c` (
  `pathway_id` int(20) NOT NULL,
  `genome_borders` varchar(20) CHARACTER SET utf8 COLLATE utf8_unicode_ci NOT NULL,
  `substrate_specificity_aa_id` int(20) NOT NULL,
  `chirality` tinyint(1) NOT NULL,
  `uni_prot_id` varchar(20) CHARACTER SET utf8 COLLATE utf8_unicode_ci NOT NULL,
  `bb_id` varchar(20) CHARACTER SET utf8 COLLATE utf8_unicode_ci NOT NULL,
  `description` text CHARACTER SET utf8 COLLATE utf8_unicode_ci NOT NULL,
  `dna_sequence` text CHARACTER SET utf8 COLLATE utf8_unicode_ci NOT NULL,
  `native_linker_before` text CHARACTER SET utf8 COLLATE utf8_unicode_ci NOT NULL,
  `native_linker_after` text CHARACTER SET utf8 COLLATE utf8_unicode_ci NOT NULL,
  PRIMARY KEY (`pathway_id`),
  KEY `substrate_specificity_aa_id` (`substrate_specificity_aa_id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8 COLLATE=utf8_unicode_ci;

--
-- Daten für Tabelle `domain_types_a_c`
--


-- --------------------------------------------------------

--
-- Tabellenstruktur für Tabelle `domain_types_e`
--

CREATE TABLE IF NOT EXISTS `domain_types_e` (
  `pathway_id` int(20) NOT NULL,
  `domain_borders` varchar(20) CHARACTER SET utf8 COLLATE utf8_unicode_ci NOT NULL,
  `uni_prot_id` varchar(20) CHARACTER SET utf8 COLLATE utf8_unicode_ci NOT NULL,
  `bb_id` varchar(20) CHARACTER SET utf8 COLLATE utf8_unicode_ci NOT NULL,
  `description` text CHARACTER SET utf8 COLLATE utf8_unicode_ci NOT NULL,
  `dna_sequence` text CHARACTER SET utf8 COLLATE utf8_unicode_ci NOT NULL,
  `native_linker_before` text CHARACTER SET utf8 COLLATE utf8_unicode_ci NOT NULL,
  `native_linker_after` text CHARACTER SET utf8 COLLATE utf8_unicode_ci NOT NULL,
  PRIMARY KEY (`pathway_id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8 COLLATE=utf8_unicode_ci;

--
-- Daten für Tabelle `domain_types_e`
--


-- --------------------------------------------------------

--
-- Tabellenstruktur für Tabelle `domain_types_others`
--

CREATE TABLE IF NOT EXISTS `domain_types_others` (
  `pathway_id` int(20) NOT NULL,
  `domain_borders` varchar(20) CHARACTER SET utf8 COLLATE utf8_unicode_ci NOT NULL,
  `type` varchar(20) CHARACTER SET utf8 COLLATE utf8_unicode_ci NOT NULL,
  `uni_prot_id` varchar(20) CHARACTER SET utf8 COLLATE utf8_unicode_ci NOT NULL,
  `bb_id` varchar(20) CHARACTER SET utf8 COLLATE utf8_unicode_ci NOT NULL,
  `description` text CHARACTER SET utf8 COLLATE utf8_unicode_ci NOT NULL,
  `dna_sequence` text CHARACTER SET utf8 COLLATE utf8_unicode_ci NOT NULL,
  `native_linker_before` text CHARACTER SET utf8 COLLATE utf8_unicode_ci NOT NULL,
  `native_linker_after` text CHARACTER SET utf8 COLLATE utf8_unicode_ci NOT NULL,
  PRIMARY KEY (`pathway_id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8 COLLATE=utf8_unicode_ci;

--
-- Daten für Tabelle `domain_types_others`
--


-- --------------------------------------------------------

--
-- Tabellenstruktur für Tabelle `domain_types_t`
--

CREATE TABLE IF NOT EXISTS `domain_types_t` (
  `pathway_id` int(20) NOT NULL,
  `domain_borders` varchar(20) CHARACTER SET utf8 COLLATE utf8_unicode_ci NOT NULL,
  `uni_prot_id` varchar(20) CHARACTER SET utf8 COLLATE utf8_unicode_ci NOT NULL,
  `bb_id` varchar(20) CHARACTER SET utf8 COLLATE utf8_unicode_ci NOT NULL,
  `description` text CHARACTER SET utf8 COLLATE utf8_unicode_ci NOT NULL,
  `dna_sequence` text CHARACTER SET utf8 COLLATE utf8_unicode_ci NOT NULL,
  `native_linker_before` text CHARACTER SET utf8 COLLATE utf8_unicode_ci NOT NULL,
  `native_linker_after` text CHARACTER SET utf8 COLLATE utf8_unicode_ci NOT NULL,
  PRIMARY KEY (`pathway_id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8 COLLATE=utf8_unicode_ci;

--
-- Daten für Tabelle `domain_types_t`
--


-- --------------------------------------------------------

--
-- Tabellenstruktur für Tabelle `domain_types_te`
--

CREATE TABLE IF NOT EXISTS `domain_types_te` (
  `pathway_id` int(20) NOT NULL,
  `domain_borders` varchar(20) CHARACTER SET utf8 COLLATE utf8_unicode_ci NOT NULL,
  `uni_prot_id` varchar(20) CHARACTER SET utf8 COLLATE utf8_unicode_ci NOT NULL,
  `bb_id` varchar(20) CHARACTER SET utf8 COLLATE utf8_unicode_ci NOT NULL,
  `description` text CHARACTER SET utf8 COLLATE utf8_unicode_ci NOT NULL,
  `dna_sequence` text CHARACTER SET utf8 COLLATE utf8_unicode_ci NOT NULL,
  `native_linker_before` text CHARACTER SET utf8 COLLATE utf8_unicode_ci NOT NULL,
  `native_linker_after` text CHARACTER SET utf8 COLLATE utf8_unicode_ci NOT NULL,
  PRIMARY KEY (`pathway_id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8 COLLATE=utf8_unicode_ci;

--
-- Daten für Tabelle `domain_types_te`
--


-- --------------------------------------------------------

--
-- Tabellenstruktur für Tabelle `main`
--

CREATE TABLE IF NOT EXISTS `main` (
  `pathway_id` int(20) NOT NULL AUTO_INCREMENT,
  `organism` varchar(20) CHARACTER SET utf8 COLLATE utf8_unicode_ci NOT NULL,
  `pathway` varchar(20) CHARACTER SET utf8 COLLATE utf8_unicode_ci NOT NULL,
  `linkout` varchar(300) CHARACTER SET utf8 COLLATE utf8_unicode_ci NOT NULL,
  `dna_sequence` text CHARACTER SET utf8 COLLATE utf8_unicode_ci NOT NULL,
  `uni_prot_id` varchar(20) CHARACTER SET utf8 COLLATE utf8_unicode_ci NOT NULL,
  `norine_id` varchar(20) CHARACTER SET utf8 COLLATE utf8_unicode_ci NOT NULL,
  PRIMARY KEY (`pathway_id`)
) ENGINE=InnoDB  DEFAULT CHARSET=utf8 COLLATE=utf8_unicode_ci AUTO_INCREMENT=2 ;

--
-- Daten für Tabelle `main`
--

INSERT INTO `main` (`pathway_id`, `organism`, `pathway`, `linkout`, `dna_sequence`, `uni_prot_id`, `norine_id`) VALUES
(1, 'Brevibacillus parabr', 'Tyrocidine', '....', '...', 'yadayada', '....');

-- --------------------------------------------------------

--
-- Tabellenstruktur für Tabelle `substrates`
--

CREATE TABLE IF NOT EXISTS `substrates` (
  `aminoacid` varchar(20) CHARACTER SET utf8 COLLATE utf8_unicode_ci NOT NULL,
  `aa_id` int(20) NOT NULL AUTO_INCREMENT,
  `structure` text CHARACTER SET utf8 COLLATE utf8_unicode_ci NOT NULL,
  `linkout` text CHARACTER SET utf8 COLLATE utf8_unicode_ci NOT NULL,
  PRIMARY KEY (`aa_id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8 COLLATE=utf8_unicode_ci AUTO_INCREMENT=1 ;

--
-- Daten für Tabelle `substrates`
--

